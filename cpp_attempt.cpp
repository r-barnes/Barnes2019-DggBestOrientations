//https://gis.stackexchange.com/questions/70800/how-to-read-vector-datasets-using-gdal-library
#include "Polygon.hpp"
#include <ogrsf_frmts.h>
#include <spatialindex/capi/sidx_api.h>
#include <spatialindex/capi/sidx_impl.h>
#include <spatialindex/capi/sidx_config.h>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <array>
#include <iomanip>
#include <cassert>
#include <fstream>
#include <chrono>

const std::string FILE_WGS84_LANDMASS = "/home/rbarnes1/scratch/dgg_best/land-polygons-complete-4326/land_polygons.shp";
const std::string FILE_OUTPUT         = "/home/rbarnes1/scratch/dgg_best/out.csv";
const std::string FILE_MERC_LANDMASS  = "/home/rbarnes1/scratch/dgg_best/land-polygons-split-3857/land_polygons.shp";

//const std::string FILE_WGS84_LANDMASS = "data/land-polygons-complete-4326/land_polygons.shp";
//const std::string FILE_OUTPUT         = "/z/out.csv";
//const std::string FILE_MERC_LANDMASS  = "data/land-polygons-split-3857/land_polygons.shp";

const double DEG_TO_RAD = M_PI/180.0;
const double RAD_TO_DEG = 180.0/M_PI;
const double IEL        = std::atan(0.5); //Icosahedron equatorial latitudes
const double IES        = 36*M_PI/180;    //Icosahedron equatorial spacing
const int PRECISION     = 1;              //Grid spacing for search

class Timer {
 private:
  typedef std::chrono::high_resolution_clock clock_;
  typedef std::chrono::duration<double, std::ratio<1> > second_;
  std::chrono::time_point<clock_> beg_;
 public:
  Timer() : beg_(clock_::now()) {}
  void reset() { beg_ = clock_::now(); }
  double elapsed() const { 
    return std::chrono::duration_cast<second_> (clock_::now() - beg_).count(); 
  }
};

template<class T>
void ToRadians(T &vec){
  for(auto &v: vec)
    v *= DEG_TO_RAD;
}

template<class T>
void ToDegrees(T &vec){
  for(auto &v: vec)
    v *= RAD_TO_DEG;
}

class LineString {
 public:
  std::vector<double> x;
  std::vector<double> y;
};

//See: https://github.com/proj4js/proj4js/blob/2006b0a06d000308caa3625005f3d5734ef11f61/lib/projections/merc.js
void WGS84toEPSG3857(
  const double lon, //Specified in radians
  const double lat, //Specified in radians
  double &x,
  double &y
){
  const double radius = 6378137; //Radius of sphere lat-lon are assumed to be on
  x = lon*radius;
  y = radius*std::log(std::tan(M_PI/4+0.5*lat));
}

template<class T>
void ToMercator(T &x, T &y){
  for(unsigned int i=0;i<x.size();i++)
    WGS84toEPSG3857(x[i],y[i],x[i],y[i]);
}

void LatLonToXYZ(
  const double lat,
  const double lon,
  const double radius,
  double &x,
  double &y,
  double &z
){
  x = radius * std::cos(lon) * std::cos(lat);
  y = radius * std::sin(lon) * std::cos(lat);
  z = radius * std::sin(lat);
}

void XYZtoLatLon(
  const double x,
  const double y,
  const double z,
  double &lat,
  double &lon
){
  double radius = x*x + y*y + z*z;
  lat = std::asin(z/radius);
  lon = std::atan2(y,x);
}

void RotatePoint(const double rlat, const double rlon, const double rtheta, double &lat, double &lon){
  lon += M_PI+rtheta;                                    //Move [0,360] and add theta
  lon  = std::fmod(2*M_PI+std::fmod(lon,2*M_PI),2*M_PI); //Map back to [0,360]
  lon -= M_PI;                                           //Move back to [-180,180] system

  double xr, yr, zr;
  LatLonToXYZ(lat, lon, 1, xr, yr, zr);
  double x =  std::cos(rlat)*std::cos(rlon)*xr + std::sin(rlon)*yr + std::sin(rlat)*std::cos(rlon)*zr;
  double y = -std::cos(rlat)*std::sin(rlon)*xr + std::cos(rlon)*yr - std::sin(rlat)*std::sin(rlon)*zr;
  double z = -std::sin(rlat)*xr + std::cos(rlat)*zr;
  XYZtoLatLon(x,y,z, lat, lon);
}


const double vla = 26.57;

class Pole {
 public:
  std::array<double,12> lon = {{     0,      0,-5*IES,-4*IES,-3*IES, -2*IES,-1*IES,   0, 1*IES,  2*IES,3*IES,4*IES}};
  std::array<double,12> lat = {{M_PI/2,-M_PI/2,   IEL,  -IEL,   IEL,   -IEL,   IEL,-IEL,   IEL,   -IEL,  IEL, -IEL}};

  Pole(){}

  void rotatePole(double rlat, double rlon, double rtheta){
    for(unsigned int i=0;i<lon.size();i++)
      RotatePoint(rlat, rlon, rtheta, lat[i], lon[i]);
  }

  void toMercator(){
    ToMercator(lon,lat);
  }

  void print() const {
    auto templon = lon;
    auto templat = lat;
    ToDegrees(templon);
    ToDegrees(templat);
    for(unsigned int i=0;i<lon.size();i++)
      std::cout<<std::setw(10)<<templat[i]<<" "<<std::setw(10)<<templon[i]<<std::endl;
  }
};

double DistancePoleToPolygon(const Pole &p, const Polygon &poly) {
  double dist = std::numeric_limits<double>::infinity();
  for(unsigned int i=0;i<p.lon.size();i++)
    dist = std::min(dist,poly.distanceFromPoint(p.lon[i],p.lat[i]));
  return dist;
}

void ReadShapefile(std::string filename, std::string layername, std::vector<Polygon> &geometries){
  GDALAllRegister();
  GDALDataset *poDS;
  poDS = (GDALDataset*) GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL);
  if( poDS == NULL ){
    std::cerr<<"Failed to open '"<<filename<<"'!"<<std::endl;
    exit( 1 );
  }

  OGRLayer *poLayer;
  poLayer = poDS->GetLayerByName(layername.c_str());
  OGRFeature *poFeature;
  poLayer->ResetReading();
  while( (poFeature = poLayer->GetNextFeature()) != NULL ){
    /*OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
    int iField;
    for( iField = 0; iField < poFDefn->GetFieldCount(); iField++ ){
      OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn( iField );
      if( poFieldDefn->GetType() == OFTInteger )
          printf( "%d,", poFeature->GetFieldAsInteger( iField ) );
      else if( poFieldDefn->GetType() == OFTInteger64 )
          printf( CPL_FRMT_GIB ",", poFeature->GetFieldAsInteger64( iField ) );
      else if( poFieldDefn->GetType() == OFTReal )
          printf( "%.3f,", poFeature->GetFieldAsDouble(iField) );
      else if( poFieldDefn->GetType() == OFTString )
          printf( "%s,", poFeature->GetFieldAsString(iField) );
      else
          printf( "%s,", poFeature->GetFieldAsString(iField) );
    }*/
    OGRGeometry *poGeometry;
    poGeometry = poFeature->GetGeometryRef();
    if (poGeometry==NULL){
      //Pass
    } else if( wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon ){
      OGRPolygon *poly = (OGRPolygon *) poGeometry;
      auto extring = poly->getExteriorRing();
      //Ignore interior rings for now: they're probably lakes
      geometries.emplace_back();
      for(int i=0;i<extring->getNumPoints();i++)
        geometries.back().exterior.emplace_back(extring->getX(i),extring->getY(i));
    } else {
      std::cerr<<"Unrecognised geometry of type: "<<wkbFlatten(poGeometry->getGeometryType())<<std::endl;
    }
    OGRFeature::DestroyFeature( poFeature );
  }
  GDALClose( poDS );
}

//Saving to disk: https://stackoverflow.com/questions/33395445/libspatialindex-loading-storing-index-on-disk
class SpIndex {
 private:
  Index* idx;
 public:
  SpIndex(){
    //Create a property set with default values.
    //See utility.cc for all defaults  http://libspatialindex.github.io/doxygen/Utility_8cc_source.html#l00031
    Tools::PropertySet* ps = GetDefaults();
    Tools::Variant var;

    //Set index type to R*-Tree
    var.m_varType   = Tools::VT_ULONG;
    var.m_val.ulVal = RT_RTree;
    ps->setProperty("IndexType", var);

    //Set index to store in memory (default is disk)
    var.m_varType   = Tools::VT_ULONG;
    var.m_val.ulVal = RT_Memory;
    ps->setProperty("IndexStorageType", var);

    //Initalise index
    idx = new Index(*ps);
    delete ps;

    //Check index is ok
    if (!idx->index().isIndexValid())
      throw std::runtime_error("Failed to create valid index");
  }

  void addBox(double left, double bottom, double right, double top, int64_t id){
    //Create array with lat/lon points
    double low[]  = {left,bottom};
    double high[] = {right,top};

    //Shapes can also have an object associated with them, but not this one!
    uint8_t* pData = 0;
    uint32_t nDataLength = 0;

    // create shape
    auto shape = SpatialIndex::Region(low,high,2);

    // insert into index along with the an object and an ID
    idx->index().insertData(nDataLength,pData,shape,id);
  }

  void addPolygon(const Polygon &p, int64_t id){
    double left   = p.minX();
    double bottom = p.minY();
    double right  = p.maxX();
    double top    = p.maxY();
    addBox(left,bottom,right,top, id);
  }

  std::vector<int64_t> overlaps(double x, double y) const {
    double coords[] = {x,y};

    //Object that will show us what results we've found
    ObjVisitor visitor;

    //Two-dimensional query point
    SpatialIndex::Point r(coords, 2);

    //Find those regions that intersect with the query point
    idx->index().intersectsWithQuery(r,visitor);

    //Copy results    
    std::vector<int64_t> temp;
    for(const auto &r: visitor.GetResults())
      temp.push_back(r->getIdentifier()); // dynamic_cast<SpatialIndex::IData*>(results[i]->clone()));

    //std::cout << "found " << visitor.GetResultCount() << " results." << std::endl;

    return temp;
  }

};


uint8_t CountOverlaps(const Pole &p, const SpIndex &sp, const std::vector<Polygon> &polygons){
  uint8_t overlaps = 0;
  for(unsigned int i=0;i<p.lat.size();i++){
    auto sp_overlaps = sp.overlaps(p.lon[i],p.lat[i]);
    if(sp_overlaps.size()==0)
      continue;
    for(auto &pi: sp_overlaps){
      if(polygons[pi].containsPoint(p.lon[i],p.lat[i])){
        overlaps++;
        break;
      }
    }
  }
  return overlaps;
}


void Test(){
  std::cerr<<"Running tests..."<<std::endl;

  {
    Pole p;
    p.rotatePole(23*DEG_TO_RAD,44*DEG_TO_RAD,5*DEG_TO_RAD);
    p.print();
  }

  SpIndex sp;
  for(double y=0;y<10;y++)
  for(double x=0;x<10;x++)
    sp.addBox(x,y,x+1,y+1,y*10+x);

  assert(sp.overlaps(3.5,3.5)[0]==33);
  assert(sp.overlaps(7.5,5.5)[0]==57);

  Polygon p;
  p.exterior.emplace_back(70,70);
  p.exterior.emplace_back(70,80);
  p.exterior.emplace_back(80,80);
  p.exterior.emplace_back(80,70);

  sp.addPolygon(p, 347);

  assert(sp.overlaps(75,77)[0]==347);

  assert(p.containsPoint(75,77));
  assert(p.containsPoint(74.3,72.2));
  assert(!p.containsPoint(69.4,72.2));

  {
    double x,y;
    WGS84toEPSG3857(-93*DEG_TO_RAD,45*DEG_TO_RAD,x,y);
    std::cout<<"(-93,45) = ("<<std::fixed<<x<<","<<std::fixed<<y<<")"<<std::endl;
    assert(std::abs(x-(-10352712.6438))<1e-4);
    assert(std::abs(y-5621521.48619)<1e-4);
  }

  std::cerr<<"Passed"<<std::endl;
}

struct POI {
  uint8_t overlaps;
  int16_t rlat;
  int16_t rlon;
  int16_t rtheta;
  double  distance;
};

std::vector<struct POI> FindPolesOfInterest(){
  std::vector<Polygon> landmass_merc;
  std::cerr<<"Reading Mercator split shapefile..."<<std::endl;
  //ReadShapefile("data/simplified-land-polygons-complete-3857/simplified_land_polygons.shp", "simplified_land_polygons", landmass_merc);
  ReadShapefile(FILE_MERC_LANDMASS, "land_polygons", landmass_merc);
  std::cerr<<"Read "<<landmass_merc.size()<<" polygons."<<std::endl;

  SpIndex sp;

  std::cerr<<"Building index..."<<std::endl;
  for(int64_t i=0;(unsigned)i<landmass_merc.size();i++) //TODO: Fix this unsigned int64 mess
    sp.addPolygon(landmass_merc[i], i);
  std::cerr<<"Index built."<<std::endl;

  std::cerr<<"Finding poles..."<<std::endl;
  std::vector<struct POI> pois;
  long count=0;

  Timer tmr;
  #pragma omp parallel for default(none) schedule(static) shared(sp,landmass_merc,pois,std::cerr) reduction(+:count)
  for(int16_t rlat  =0; rlat  <634; rlat  +=PRECISION) //(pi()/2-vla)*180/pi()
  for(int16_t rlon  =0; rlon  <720; rlon  +=PRECISION)
  for(int16_t rtheta=0; rtheta<720; rtheta+=PRECISION){
    count++;
    Pole p;
    p.rotatePole(rlat/10.0*DEG_TO_RAD, rlon/10.0*DEG_TO_RAD, rtheta/10.0*DEG_TO_RAD);
    p.toMercator();
    auto overlaps = CountOverlaps(p,sp,landmass_merc);
    if(overlaps==0 || overlaps>=8){
      #pragma omp critical
      pois.push_back(POI{overlaps,rlat,rlon,rtheta,-1});
    }
  }

  double t = tmr.elapsed();
  std::cout << "Time taken = " << t <<"s"<< std::endl;

  std::cerr<<"Checked = "<<count<<std::endl;
  std::cerr<<"Found "<<pois.size()<<" poles of interest."<<std::endl;

  return pois;
}

void DistancesToPoles(std::vector<struct POI> &pois){
  std::cerr<<"Reading WGS84 shapefile..."<<std::endl;
  std::vector<Polygon> landmass_wgs84;
  ReadShapefile(FILE_WGS84_LANDMASS, "land_polygons", landmass_wgs84);

  for(auto &p: landmass_wgs84)
    p.toRadians();

  Timer tmr;

  std::cerr<<"Calculating distances to poles..."<<std::endl;
  #pragma omp parallel for default(none) shared(pois,landmass_wgs84)
  for(unsigned int i=0;i<pois.size();i++){
    Pole p;
    p.rotatePole(pois[i].rlat/10.0*DEG_TO_RAD, pois[i].rlon/10.0*DEG_TO_RAD, pois[i].rtheta/10.0*DEG_TO_RAD);
    pois[i].distance = std::numeric_limits<double>::infinity();
    for(const auto &g: landmass_wgs84)
      pois[i].distance = std::min(pois[i].distance, DistancePoleToPolygon(p,g));
  }

  double t = tmr.elapsed();
  std::cout << "Time taken = " << t <<"s"<< std::endl;
}

int main(int argc, char **argv){
  Test();

  std::cerr<<"PRECISION = "<<PRECISION<<std::endl;

  auto pois = FindPolesOfInterest();

  DistancesToPoles(pois);

  std::cerr<<"Sorting poles of interest..."<<std::endl;
  std::sort(pois.begin(),pois.end(), [](const POI &a, const POI &b){
    return a.distance>b.distance;
  });

  std::cerr<<"Writing output..."<<std::endl;
  std::ofstream fout(FILE_OUTPUT);
  for(const auto &p: pois)
    fout<<(int)p.overlaps<<","<<p.rlat<<","<<p.rlon<<","<<p.rtheta<<","<<p.distance<<"\n";
  fout.close();
}
