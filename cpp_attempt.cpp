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
  lon += M_PI+rtheta;                          //Move [0,360] and add theta
  lon  = std::fmod(2*M_PI+std::fmod(lon,2*M_PI),2*M_PI); //Map back to [0,360]
  lon -= M_PI;                                 //Move back to [-180,180] system

  double xr, yr, zr;
  LatLonToXYZ(lat, lon, 1, xr, yr, zr);
  double x =  std::cos(rlat)*std::cos(rlon)*xr + std::sin(rlon)*yr + std::sin(rlat)*std::cos(rlon)*zr;
  double y = -std::cos(rlat)*std::sin(rlon)*xr + std::cos(rlon)*yr - std::sin(rlat)*std::sin(rlon)*zr;
  double z = -std::sin(rlat)*xr + std::cos(rlat)*zr;
  XYZtoLatLon(x,y,z, lat, lon);
}

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
};



const double vla = 26.57;

class Pole {
 public:
  std::array<double,12> plon = {{ 90,0,-180,-144,-108, -72,-36,   0, 36,  72,108, 144}};
  std::array<double,12> plat = {{-90,0, vla,-vla, vla,-vla,vla,-vla,vla,-vla,vla,-vla}};
  Pole(){
    ToRadians(plon);
    ToRadians(plat);
  }

  void rotatePole(double rlat, double rlon, double rtheta){
    rlat *= DEG_TO_RAD;
    rlon *= DEG_TO_RAD;

    for(int i=0;i<plon.size();i++)
      RotatePoint(rlat, rlon, rtheta, plat[i], plon[i]);
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
    double low[]    = {left,bottom};
    double high[]   = {right,top};

    //Shapes can also have an object associated with them, but not this one!
    uint8_t* pData = 0;
    uint32_t nDataLength = 0;

    // create shape
    SpatialIndex::IShape* shape = new SpatialIndex::Region(low,high,2);

    // insert into index along with the an object and an ID
    idx->index().insertData(nDataLength,pData,*shape,id);

    delete shape;

  }

  void addPolygon(const Polygon &p, int64_t id){
    double left   = p.minX();
    double bottom = p.minY();
    double right  = p.maxX();
    double top    = p.maxY();
    addBox(left,bottom,right,top, id);
  }

  std::vector<int64_t> overlaps(double x, double y, double maxResults){
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





int main(int argc, char **argv){
  std::vector<Polygon> geometries;
  ReadShapefile("data/simplified-land-polygons-complete-3857/simplified_land_polygons.shp", "simplified_land_polygons", geometries);
  //ReadShapefile("data/land-polygons-split-3857/land_polygons.shp", "land_polygons", geometries);
  std::cerr<<"Read "<<geometries.size()<<" polygons."<<std::endl;

  SpIndex sp;

  std::cerr<<"Building index..."<<std::endl;
  for(int64_t i=0;i<geometries.size();i++)
    sp.addPolygon(geometries[i], i);
  std::cerr<<"Index built."<<std::endl;

  for(auto &i: sp.overlaps(-11501094,3771603,5000))
    std::cerr<<"Found "<<i<<std::endl;

  Pole p;
  p.rotatePole(23*DEG_TO_RAD,44*DEG_TO_RAD,5*DEG_TO_RAD);
}






















/*
  void toMercator(){
    // EPSG:4326 definition: http://spatialreference.org/ref/epsg/4326/proj4/
    const char* wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";
    // EPSG:27700 definition: http://spatialreference.org/ref/epsg/27700/proj4/
    const char* merc = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs";
    
    projPJ pj_src = pj_init_plus(src);
    projPJ pj_dst = pj_init_plus(dst);

    toRadian();

    pj_transform(wgs84, merc, exterior.x.size(), 1, exterior.x.data(), exterior.y.data(), NULL);
    pj_transform(wgs84, merc, interior.x.size(), 1, interior.x.data(), interior.y.data(), NULL);
  }
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