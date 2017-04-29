//https://gis.stackexchange.com/questions/70800/how-to-read-vector-datasets-using-gdal-library
#include <gdal/ogrsf_frmts.h>
#include <proj_api.h>
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

const double IEL    = std::atan(0.5); //Icosahedron equatorial latitudes
const double IES    = 36*M_PI/180;    //Icosahedron equatorial spacing
const int PRECISION = 50;             //Grid spacing for search

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

const char* wgs84_str = "+init=epsg:4326"; // EPSG:4326 definition: http://spatialreference.org/ref/epsg/4326/proj4/
const char* merc_str  = "+init=epsg:3857";

//Using static means that these variables are file-scope. If we didn't use
//static then the variables might be accessed across many files and OpenMP would
//not like that
static projPJ  pj_wgs84;
static projPJ  pj_merc;
static projCtx pj_ctx;

//Set it so each thread has its own, private copy of these variables
#pragma omp threadprivate(pj_wgs84)
#pragma omp threadprivate(pj_merc)
#pragma omp threadprivate(pj_ctx)


template<class T>
void ToMercator(T &x, T &y){
  pj_transform(pj_wgs84, pj_merc, x.size(), 1, x.data(), y.data(), NULL);
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

//Info: https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
//Info: https://stackoverflow.com/questions/8721406/how-to-determine-if-a-point-is-inside-a-2d-convex-polygon/23223947#23223947
//Info: http://stackoverflow.com/a/2922778/752843
int ContainsPoint(const std::vector<double> &vertx, const std::vector<double> &verty, const double testx, const double testy) {
  unsigned int i, j;
  int c = 0;
  for (i = 0, j = vertx.size()-1; i < vertx.size(); j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
     (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}

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
  const double Rearth = 6371; //km

  return Rearth*std::acos( std::sin(lat1)*std::sin(lat2) + std::cos(lat1)*std::cos(lat2)*std::cos(lon2-lon1) );
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


class Polygon {
 public:
  LineString interior;
  LineString exterior;
  double minX() const {
    return *std::min_element(exterior.x.begin(), exterior.x.end());
  }
  double maxX() const {
    return *std::max_element(exterior.x.begin(), exterior.x.end());
  }
  double minY() const {
    return *std::min_element(exterior.y.begin(), exterior.y.end());
  }
  double maxY() const {
    return *std::max_element(exterior.y.begin(), exterior.y.end());
  }
  void toRadian() {
    for(auto &i: exterior.x) i *= DEG_TO_RAD;
    for(auto &i: exterior.y) i *= DEG_TO_RAD;
    for(auto &i: interior.x) i *= DEG_TO_RAD;
    for(auto &i: interior.y) i *= DEG_TO_RAD;
  }
  void toDegree() {
    for(auto &i: exterior.x) i *= RAD_TO_DEG;
    for(auto &i: exterior.y) i *= RAD_TO_DEG;
    for(auto &i: interior.x) i *= RAD_TO_DEG;
    for(auto &i: interior.y) i *= RAD_TO_DEG;
  }
  bool containsPoint(const double testx, const double testy) const {
    return ContainsPoint(exterior.x, exterior.y, testx, testy);
  }

  // """
  // Calculate the closest distance between a polygon and a latitude-longitude
  // point, using only spherical considerations. Ignore edges.
  // :param lat  Latitude of query point in degrees
  // :param lon  Longitude of query point in degrees
  // :param geom A `shapely` geometry whose points are in latitude-longitude space
  // :returns: The minimum distance in kilometres between the polygon and the
  //           query point
  // """
  double distanceFromPoint(const double px, const double py) const {
    double dist = std::numeric_limits<double>::infinity();
    for(unsigned int i=0;i<exterior.x.size();i++)
      dist = std::min(dist,GeoDistance(px,py,exterior.x[i],exterior.y[i]));
    return dist;
  }

  double distanceFromPole(const Pole &p) const {
    double dist = std::numeric_limits<double>::infinity();
    for(unsigned int i=0;i<p.lon.size();i++)
      dist = std::min(dist,distanceFromPoint(p.lon[i],p.lat[i]));
    return dist;
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
      for(int i=0;i<extring->getNumPoints();i++){
        geometries.back().exterior.x.emplace_back(extring->getX(i));
        geometries.back().exterior.y.emplace_back(extring->getY(i));
      }
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

  std::cerr<<"PROJ version = "<<PJ_VERSION<<std::endl;

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
  p.exterior.x.push_back(70); p.exterior.y.push_back(70);
  p.exterior.x.push_back(70); p.exterior.y.push_back(80);
  p.exterior.x.push_back(80); p.exterior.y.push_back(80);
  p.exterior.x.push_back(80); p.exterior.y.push_back(70);

  sp.addPolygon(p, 347);

  assert(sp.overlaps(75,77)[0]==347);

  assert(p.containsPoint(75,77));
  assert(p.containsPoint(74.3,72.2));
  assert(!p.containsPoint(69.4,72.2));

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
  ReadShapefile("data/simplified-land-polygons-complete-3857/simplified_land_polygons.shp", "simplified_land_polygons", landmass_merc);
  //ReadShapefile("data/land-polygons-split-3857/land_polygons.shp", "land_polygons", landmass_merc);
  std::cerr<<"Read "<<landmass_merc.size()<<" polygons."<<std::endl;

  SpIndex sp;

  std::cerr<<"Building index..."<<std::endl;
  for(int64_t i=0;(unsigned)i<landmass_merc.size();i++) //TODO: Fix this unsigned int64 mess
    sp.addPolygon(landmass_merc[i], i);
  std::cerr<<"Index built."<<std::endl;

  std::cerr<<"Finding poles..."<<std::endl;
  std::vector<struct POI> pois;
  #pragma omp parallel for collapse(3)
  for(int16_t rlat  =0; rlat  <634; rlat  +=PRECISION) //(pi()/2-vla)*180/pi()
  for(int16_t rlon  =0; rlon  <720; rlon  +=PRECISION)
  for(int16_t rtheta=0; rtheta<720; rtheta+=PRECISION){
    Pole p;
    p.rotatePole(rlat/10.0*DEG_TO_RAD, rlon/10.0*DEG_TO_RAD, rtheta/10.0*DEG_TO_RAD);
    p.toMercator();
    auto overlaps = CountOverlaps(p,sp,landmass_merc);
    if(overlaps==0 || overlaps>=8){
      #pragma omp critical
      pois.push_back(POI{overlaps,rlat,rlon,rtheta,-1});
    }
  }

  std::cerr<<"Found "<<pois.size()<<" poles of interest."<<std::endl;

  return pois;
}

void DistancesToPoles(std::vector<struct POI> &pois){
  std::cerr<<"Reading WGS84 shapefile..."<<std::endl;
  std::vector<Polygon> landmass_wgs84;
  ReadShapefile("data/land-polygons-complete-4326/land_polygons.shp", "land_polygons", landmass_wgs84);

  std::cerr<<"Calculating distances to poles..."<<std::endl;
  for(auto &poi: pois){
    Pole p;
    p.rotatePole(poi.rlat*DEG_TO_RAD, poi.rlon*DEG_TO_RAD, poi.rtheta*DEG_TO_RAD);
    poi.distance = std::numeric_limits<double>::infinity();
    for(const auto &g: landmass_wgs84)
      poi.distance = std::min(poi.distance, g.distanceFromPole(p));
  }
}

int main(int argc, char **argv){
  Test();

  auto pois = FindPolesOfInterest();

  //DistancesToPoles(pois);
}























/*

  void toWGS84(){
    // EPSG:4326 definition: http://spatialreference.org/ref/epsg/4326/proj4/
    const char* wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";
    // EPSG:27700 definition: http://spatialreference.org/ref/epsg/27700/proj4/
    const char* merc = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs";
    
    projPJ pj_src = pj_init_plus(src);
    projPJ pj_dst = pj_init_plus(dst);

    toRadian();

    pj_transform(merc, wgs84, exterior.x.size(), 1, exterior.x.data(), exterior.y.data(), NULL);
    pj_transform(merc, wgs84, interior.x.size(), 1, interior.x.data(), interior.y.data(), NULL);
  }
  */