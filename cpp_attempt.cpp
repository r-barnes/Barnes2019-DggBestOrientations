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
#include <bitset>

#ifdef ENV_XSEDE
  const std::string FILE_WGS84_LANDMASS = "/home/rbarnes1/scratch/dgg_best/land-polygons-complete-4326/land_polygons.shp";
  const std::string FILE_OUTPUT_ROT     = "/home/rbarnes1/scratch/dgg_best/out-rot.csv";
  const std::string FILE_OUTPUT_VERT    = "/home/rbarnes1/scratch/dgg_best/out-vert.csv";
  const std::string FILE_MERC_LANDMASS  = "/home/rbarnes1/scratch/dgg_best/land-polygons-split-3857/land_polygons.shp";
#elif ENV_LAPTOP
  const std::string FILE_WGS84_LANDMASS = "data/land-polygons-complete-4326/land_polygons.shp";
  const std::string FILE_OUTPUT_ROT     = "/z/out-rot.csv";
  const std::string FILE_OUTPUT_VERT    = "/z/out-vert.csv";
  const std::string FILE_MERC_LANDMASS  = "data/land-polygons-split-3857/land_polygons.shp";
#else
  this-is-an-error
#endif

const double DEG_TO_RAD = M_PI/180.0;
const double RAD_TO_DEG = 180.0/M_PI;
const double IEL        = std::atan(0.5); //Icosahedron equatorial latitudes
const double IES        = 36*M_PI/180;    //Icosahedron equatorial spacing
const double PRECISION  = 0.05;           //Grid spacing for search
const double DIV        = 100.0;

// """
// Calculate the Great Circle distance on Earth between two latitude-longitude
// points
// :param lat1 Latitude of Point 1 in degrees
// :param lon1 Longtiude of Point 1 in degrees
// :param lat2 Latitude of Point 2 in degrees
// :param lon2 Longtiude of Point 2 in degrees
// :returns Distance between the two points in kilometres
// """
double GeoDistanceFlatEarth(
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
}

double GeoDistanceHaversine(
  const double lon1, 
  const double lat1, 
  const double lon2,
  const double lat2 
){
  const double Rearth = 6371; //km
  const double dlon   = lon2 - lon1;
  const double dlat   = lat2 - lat1;
  const double a      = std::pow(std::sin(dlat/2),2) + std::cos(lat1) * std::cos(lat2) * std::pow(std::sin(dlon/2),2);
  const double c      = 2*std::asin(std::sqrt(a));
  return Rearth*c;
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


class Pole {
 public:
  std::array<double,12> lon = {{     0,      0,-5*IES,-4*IES,-3*IES, -2*IES,-1*IES,   0, 1*IES,  2*IES,3*IES,4*IES}};
  std::array<double,12> lat = {{M_PI/2,-M_PI/2,   IEL,  -IEL,   IEL,   -IEL,   IEL,-IEL,   IEL,   -IEL,  IEL, -IEL}};

  Pole(){}

  Pole(double rlat, double rlon, double rtheta){
    rotatePole(rlat,rlon,rtheta);
  }

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

void ReadShapefile(std::string filename, std::string layername, Polygons &geometries){
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

std::bitset<12> CountOverlaps(const Pole &p, const SpIndex &sp, const std::vector<Polygon> &polygons){
  std::bitset<12> overlaps = 0;
  for(unsigned int i=0;i<p.lat.size();i++){
    if(p.lat[i]>83.7*DEG_TO_RAD) //The island "83-42" as at 83.7N so anything north of this is on water
      continue;
    if(p.lat[i]<-80*DEG_TO_RAD){ //The southmost extent of water is ~79.5S, so anything south of this is on land
      overlaps.set(i);
      continue;
    }
    double x,y;
    WGS84toEPSG3857(p.lon[i],p.lat[i],x,y);
    const auto pid = sp.queryPoint(x,y);
    if(pid==-1)
      continue;
    if(polygons.at(pid).containsPoint(x,y))
      overlaps.set(i);
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
    p.rotatePole(22.8*DEG_TO_RAD,3.6*DEG_TO_RAD,45.6*DEG_TO_RAD);
    p.print();
    p.toMercator();
    p.print();
  }

  SpIndex sp;
  int id=0;
  for(double y=0;y<1000;y+=100)
  for(double x=0;x<1000;x+=100)
    sp.addBox(x,y,x+100,y+100,id++);

  assert(sp.queryPoint(350,350)==33);
  assert(sp.queryPoint(750,550)==57);

  Polygon p;
  p.exterior.emplace_back(1200,1200);
  p.exterior.emplace_back(1200,1300);
  p.exterior.emplace_back(1300,1300);
  p.exterior.emplace_back(1300,1200);

  AddPolygonToSpIndex(p, sp, 347);
  sp.buildIndex();

  std::cerr<<sp.queryPoint(1250,1270)<<std::endl;

  assert(sp.queryPoint(1250,1270)==347);

  assert(p.containsPoint(1250,1270));
  assert(p.containsPoint(1243,1222));
  assert(!p.containsPoint(1194,1222));

  {
    double x,y;
    WGS84toEPSG3857(-93*DEG_TO_RAD,45*DEG_TO_RAD,x,y);
    std::cout<<"(-93,45) = ("<<std::fixed<<x<<","<<std::fixed<<y<<")"<<std::endl;
    assert(std::abs(x-(-10352712.6438))<1e-4);
    assert(std::abs(y-5621521.48619)<1e-4);
  }

  std::cerr<<"Passed"<<std::endl;
}

class POI {
 public:
  std::bitset<12> overlaps = 0;
  int16_t         rlat     = 0;
  int16_t         rlon     = 0;
  int16_t         rtheta   = 0;
  double          mindist  = std::numeric_limits<double>::infinity();
  double          maxdist  = -std::numeric_limits<double>::infinity();
  double          avgdist  = 0;
  uint16_t        cluster  = 0;
  POI(std::bitset<12> overlaps,double rlat,double rlon,double rtheta){
    this->overlaps = overlaps;
    this->rlat     = rlat;
    this->rlon     = rlon;
    this->rtheta   = rtheta;
  }
};

std::vector<struct POI> FindPolesOfInterest(
  const Polygons &landmass_merc,
  const SpIndex &sp
){
  std::cerr<<"Finding poles..."<<std::endl;
  std::vector<struct POI> pois;
  long count=0;

  Timer tmr;
  #pragma omp parallel for default(none) schedule(static) shared(sp,landmass_merc,pois,std::cerr) reduction(+:count)
  for(int16_t rlat  =0; rlat  <(int)(63.4*DIV); rlat  +=(int)(PRECISION*DIV)) //(pi()/2-IEL)*180/pi()
  for(int16_t rlon  =0; rlon  <(int)(72.0*DIV); rlon  +=(int)(PRECISION*DIV))
  for(int16_t rtheta=0; rtheta<(int)(72.0*DIV); rtheta+=(int)(PRECISION*DIV)){
    count++;
    Pole p(rlat/DIV*DEG_TO_RAD, rlon/DIV*DEG_TO_RAD, rtheta/DIV*DEG_TO_RAD);
    auto overlaps = CountOverlaps(p,sp,landmass_merc);  
    if(overlaps==0 || overlaps.count()>=8){
      #pragma omp critical
      pois.emplace_back(overlaps,rlat,rlon,rtheta);
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
  for(unsigned int pn=0;pn<pois.size();pn++){
    Pole p;
    p.rotatePole(pois[pn].rlat/DIV*DEG_TO_RAD, pois[pn].rlon/DIV*DEG_TO_RAD, pois[pn].rtheta/DIV*DEG_TO_RAD);

    for(unsigned int j=0;j<p.lat.size();j++){
      double xr,yr,zr;
      LatLonToXYZ(p.lat[j],p.lon[j],1,xr,yr,zr);
      const auto cp = pc.queryPoint(xr,yr,zr); //Closest point
      double plat,plon;
      XYZtoLatLon(cp.x,cp.y,cp.z,plat,plon);
      auto dist = GeoDistanceFlatEarth(plon,plat,p.lon[j],p.lat[j]);
      if(pois[pn].overlaps.test(j))
        dist = -dist;
      pois[pn].mindist = std::min(pois[pn].mindist,dist);
      pois[pn].maxdist = std::max(pois[pn].maxdist,dist);
      pois[pn].avgdist += dist;
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
    return a.mindist>b.mindist;
  });

  std::cerr<<"Writing output..."<<std::endl;
  {
    std::ofstream fout(FILE_OUTPUT_ROT);
    fout<<"Num,Cluster,Overlaps,OverlapCount,Lat,Lon,Theta,MinDistance,MaxDistance,AvgDistance\n";
    for(unsigned int i=0;i<pois.size();i++)
      fout<<i<<","
          <<pois[i].cluster             <<","
          <<pois[i].overlaps.to_string()<<","
          <<pois[i].overlaps.count()    <<","
          <<pois[i].rlat                <<","
          <<pois[i].rlon                <<","
          <<pois[i].rtheta              <<","
          <<pois[i].mindist             <<","
          <<pois[i].maxdist             <<","
          <<(pois[i].avgdist/12)
          <<"\n";
  }

  {
    std::ofstream fout(FILE_OUTPUT_VERT);
    fout<<"Num,Cluster,Lat,Lon,OnLand\n";
    for(unsigned int pn=0;pn<pois.size();pn++){
      Pole pole(pois[pn].rlat/DIV*DEG_TO_RAD,pois[pn].rlon/DIV*DEG_TO_RAD,pois[pn].rtheta/DIV*DEG_TO_RAD);
      for(unsigned int i=0;i<pole.lat.size();i++)
        fout<<pn                    <<","
            <<pois[pn].cluster      <<","
            <<pole.lat[i]*RAD_TO_DEG<<","
            <<pole.lon[i]*RAD_TO_DEG<<","
            <<pois[pn].overlaps.test(i)
            <<"\n";
    }
  }


}
