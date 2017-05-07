//https://gis.stackexchange.com/questions/70800/how-to-read-vector-datasets-using-gdal-library
//https://math.stackexchange.com/questions/44080/approximating-geodesic-distances-on-the-sphere-by-euclidean-distances-of-a-trans
#include "Polygon.hpp"
#include "SpIndex.hpp"
#include "PointCloud.hpp"
#include <ogrsf_frmts.h>
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

// """
// Calculate the Great Circle distance on Earth between two latitude-longitude
// points
// :param lat1 Latitude of Point 1 in degrees
// :param lon1 Longtiude of Point 1 in degrees
// :param lat2 Latitude of Point 2 in degrees
// :param lon2 Longtiude of Point 2 in degrees
// :returns Distance between the two points in kilometres
// """
double GeoDistance(
  const double lon1, 
  const double lat1, 
  const double lon2,
  const double lat2 
){
  //Flat Earth Approx
  const double Rearth = 6371; //km
  const double dlat = lat2-lat1;
  const double dlon = lon2-lon1;
  const double mlat = (lat2+lat1)/2;
  return Rearth*std::sqrt(dlat*dlat + std::pow(std::cos(mlat)*dlon,2));

  //https://www.mapbox.com/blog/cheap-ruler/
  //http://www.focusonmath.org/sites/focusonmath.org/files/assets/MT2004-08-20a%281%29.pdf
  //From: https://www.gpo.gov/fdsys/pkg/CFR-2005-title47-vol4/pdf/CFR-2005-title47-vol4-sec73-208.pdf. Valid for distances <295 miles
  // const double ml      = (lat1+lat2)/2;
  // const double cs      = std::cos(ml);
  // const double cs2     = 2*cs*cs  - 1;
  // const double cs3     = 2*cs*cs2 - cs;
  // const double cs4     = 2*cs*cs3 - cs2;
  // const double cs5     = 2*cs*cs4 - cs3;
  // const double kpd_lat = 111.13209    - 0.56605*cs2 + 0.00120*cs4;
  // const double kpd_lon = 111.41513*cs - 0.09455*cs3 + 0.00012*cs5;
  // const double ns      = kpd_lat*(lat1-lat2);
  // const double ew      = kpd_lon*(lon1-lon2);
  // return std::sqrt(ns*ns+ew*ew);

  //Law of Cosines method
  //const double Rearth = 6371; //km
  //return Rearth*std::acos( std::sin(lat1)*std::sin(lat2) + std::cos(lat1)*std::cos(lat2)*std::cos(lon2-lon1) );

  //Haversine Distance
  // const double Rearth = 6371; //km
  // const double dlon   = lon2 - lon1;
  // const double dlat   = lat2 - lat1;
  // const double a      = std::pow(std::sin(dlat/2),2) + std::cos(lat1) * std::cos(lat2) * std::pow(std::sin(dlon/2),2);
  // const double c      = 2*std::asin(std::sqrt(a));
  // return Rearth*c;
}

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

//See: https://github.com/proj4js/proj4js/blob/2006b0a06d000308caa3625005f3d5734ef11f61/lib/projections/merc.js
void WGS84toEPSG3857(
  const double lon, //Specified in radians
  const double lat, //Specified in radians
  double &x,
  double &y
){
  //Radius of sphere lat-lon are assumed to be on (radius of Earth in metres).
  //This has been chosen to match the OpenStreetMap Mercator projection. Change
  //with care or point-in-polygon testing may fail.
  const double radius = 6378137; 
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
    for(unsigned int i=0;i<lon.size();i++)
      std::cout<<std::fixed<<std::setw(10)<<templat[i]<<" "<<std::fixed<<std::setw(10)<<templon[i]<<" -- "<<std::setw(10)<<std::fixed<<templat[i]*RAD_TO_DEG<<" "<<std::setw(10)<<std::fixed<<templon[i]*RAD_TO_DEG<<std::endl;
  }
};

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

uint8_t CountOverlaps(const Pole &p, SpIndex &sp, const std::vector<Polygon> &polygons){
  uint8_t overlaps = 0;
  for(unsigned int i=0;i<p.lat.size();i++){
    const auto pid = sp.queryPoint(p.lon[i],p.lat[i]);
    if(pid==-1)
      continue;
    if(polygons.at(pid).containsPoint(p.lon[i],p.lat[i]))
      overlaps++;
  }
  return overlaps;
}

void AddPolygonToSpIndex(const Polygon &poly, SpIndex &sp, const int id){
  const int xmin = poly.minX();
  const int ymin = poly.minY();
  const int xmax = poly.maxX();
  const int ymax = poly.maxY();
  sp.addBoxDeferred(xmin,ymin,xmax,ymax,id);
}

void Test(){
  std::cerr<<"Running tests..."<<std::endl;

  {
    Pole p;
    p.rotatePole(23*DEG_TO_RAD,44*DEG_TO_RAD,5*DEG_TO_RAD);
    p.print();
    p.toMercator();
    p.print();
  }

  SpIndex sp;
  int id=0;
  for(double y=0;y<100;y+=10)
  for(double x=0;x<100;x+=10)
    sp.addBox(x,y,x+10,y+10,id++);

  assert(sp.queryPoint(35,35)==33);
  assert(sp.queryPoint(75,55)==57);

  Polygon p;
  p.exterior.emplace_back(120,120);
  p.exterior.emplace_back(120,130);
  p.exterior.emplace_back(130,130);
  p.exterior.emplace_back(130,120);

  AddPolygonToSpIndex(p, sp, 347);

  assert(sp.queryPoint(125,127)==347);

  assert(p.containsPoint(125,127));
  assert(p.containsPoint(124.3,122.2));
  assert(!p.containsPoint(119.4,122.2));

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

  {
    std::cerr<<"Building index..."<<std::endl;
    Timer tmr;
    for(int64_t i=0;(unsigned)i<landmass_merc.size();i++)
      AddPolygonToSpIndex(landmass_merc[i], sp, i);
    sp.buildIndex();
    std::cerr<<"Index built."<<std::endl;
    std::cerr<<"Time taken = "<<tmr.elapsed()<<std::endl;
  }

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

  PointCloud pc;

  for(auto &p: landmass_wgs84)
    p.toRadians();

  std::cerr<<"Adding polygon points to PointCloud..."<<std::endl;
  for(const auto &poly: landmass_wgs84){
    for(const auto &p: poly.exterior){
      double xr, yr, zr;
      LatLonToXYZ(p.y, p.x, 1, xr, yr, zr);
      pc.addPoint(xr,yr,zr);
    }
  }

  std::cerr<<"Building kd-tree..."<<std::endl;
  pc.buildIndex();

  Timer tmr;

  std::cerr<<"Calculating distances to poles..."<<std::endl;
  #pragma omp parallel for default(none) shared(pois,pc)
  for(unsigned int i=0;i<pois.size();i++){
    Pole p;
    p.rotatePole(pois[i].rlat/10.0*DEG_TO_RAD, pois[i].rlon/10.0*DEG_TO_RAD, pois[i].rtheta/10.0*DEG_TO_RAD);

    pois[i].distance = std::numeric_limits<double>::infinity();
    for(unsigned int j=0;j<p.lat.size();j++){
      double xr,yr,zr;
      LatLonToXYZ(p.lat[j],p.lon[j],1,xr,yr,zr);
      const auto cp = pc.queryPoint(xr,yr,zr); //Closest point
      double plat,plon;
      XYZtoLatLon(cp.x,cp.y,cp.z,plat,plon);
      const auto dist = GeoDistance(plon,plat,p.lon[j],p.lat[j]);
      pois[i].distance = std::min(pois[i].distance,dist);
    }
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
