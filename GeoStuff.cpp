#include "GeoStuff.hpp"
#include <cmath>
#include "doctest.h"
#include <ogrsf_frmts.h>
#include <cmath>
#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/Constants.hpp>

const double RAD_TO_DEG = 180.0/M_PI;

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
  const Point2D &a,
  const Point2D &b
){
  //Flat Earth Approx
  const double Rearth = 6371; //km
  const double dlat = b.y-a.y;
  const double dlon = b.x-a.x;
  const double mlat = (b.y+a.y)/2;
  return Rearth*std::sqrt(dlat*dlat + std::pow(std::cos(mlat)*dlon,2));
}

TEST_CASE("GeoDistanceFlatEarth"){
  auto sanfran  = Point2D(-122.4194, 37.7749).toRadians();
  auto berkeley = Point2D(-122.2727, 37.8716).toRadians();
  auto dist = GeoDistanceFlatEarth(sanfran, berkeley);
  CHECK(dist==doctest::Approx(16.782).epsilon(0.01));
}



double GeoDistanceHaversine(
  const Point2D &pa,
  const Point2D &pb
){
  const double Rearth = 6371; //km
  const double dlon   = pb.x-pa.x;
  const double dlat   = pb.y-pa.y;
  const double a      = std::pow(std::sin(dlat/2),2) + std::cos(pa.y) * std::cos(pb.y) * std::pow(std::sin(dlon/2),2);
  const double c      = 2*std::asin(std::sqrt(a));
  return Rearth*c;
}

TEST_CASE("GeoDistanceHaversine"){
  auto seattle = Point2D(-122.3321, 46.6062).toRadians();
  auto johan   = Point2D(28.0473, 26.2041).toRadians();
  auto dist = GeoDistanceHaversine(seattle,johan);
  CHECK(dist==doctest::Approx(11387.97));
}



double GeoDistanceEllipsoid(
  const Point2D &a,
  const Point2D &b
){
  static auto geod = GeographicLib::Geodesic(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f());

  const auto gline = geod.InverseLine(
    a.y*RAD_TO_DEG,
    a.x*RAD_TO_DEG,
    b.y*RAD_TO_DEG,
    b.x*RAD_TO_DEG
  );

  return gline.Distance()/1000; //km
}



//See: https://github.com/proj4js/proj4js/blob/2006b0a06d000308caa3625005f3d5734ef11f61/lib/projections/merc.js
Point2D WGS84toEPSG3857(
  const Point2D &p //Specified in radians
){
  //Radius of sphere lat-lon are assumed to be on (radius of Earth in metres).
  //This has been chosen to match the OpenStreetMap Mercator projection. Change
  //with care or point-in-polygon testing may fail.
  const double radius = 6378137;
  return Point2D (
    p.x*radius,
    radius*std::log(std::tan(M_PI/4+0.5*p.y))
  );
}

TEST_CASE("WGS84toEPSG3857"){
  SUBCASE("Minneapolis"){
    const auto a = WGS84toEPSG3857(Point2D(-93.2167,44.8803).toRadians());
    CHECK(a.x==doctest::Approx(-10376835.57891));
    CHECK(a.y==doctest::Approx(5602696.81363));
  }
  SUBCASE("Fairbanks"){
    const auto a = WGS84toEPSG3857(Point2D(-147.856,64.815).toRadians());
    CHECK(a.x==doctest::Approx(-16459254.63074096));
    CHECK(a.y==doctest::Approx(9559809.55561878));
  }
  SUBCASE("Chatham Island, New Zealand"){
    const auto a = WGS84toEPSG3857(Point2D(-176.457,-43.81).toRadians());
    CHECK(a.x==doctest::Approx(-19643103));
    CHECK(a.y==doctest::Approx(-5436086));
  }
  SUBCASE("Anadyr, Russia"){
    const auto a = WGS84toEPSG3857(Point2D(177.741,64.7347).toRadians());
    CHECK(a.x==doctest::Approx(19786038));
    CHECK(a.y==doctest::Approx(9538835));
  }
}



double EuclideanDistance(const Point3D &a, const Point3D &b){
  return std::sqrt(std::pow(a.x-b.x,2)+std::pow(a.y-b.y,2)+std::pow(a.z-b.z,2));
}

TEST_CASE("EuclideanDistance"){
  auto a = Point3D(10, 33, 17);
  auto b = Point3D(97, 42, 3);
  auto dist = EuclideanDistance(a,b);
  CHECK(dist == doctest::Approx(88.5776495));
}



Point3D WGS84toEllipsoidCartesian(const Point2D &p) {
  static const GeographicLib::Geocentric& earth = GeographicLib::Geocentric::WGS84();

  Point3D temp;

  //earth.Forward(lat, lon, h, X, Y, Z);
  earth.Forward(p.y*RAD_TO_DEG,p.x*RAD_TO_DEG,0,temp.x,temp.y,temp.z);

  return temp;
}

Point2D EllipsoidCartesiantoWGS84(const Point3D &p){
  static const GeographicLib::Geocentric& earth = GeographicLib::Geocentric::WGS84();

  Point2D temp;
  double height;
  earth.Reverse(p.x,p.y,p.z,temp.y,temp.x,height);
  CHECK(height==doctest::Approx(0));
  temp.toRadians();
  return temp;
}



Polygons ReadShapefile(std::string filename, std::string layername){
  GDALAllRegister();
  GDALDataset *poDS;
  poDS = (GDALDataset*) GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL);
  if( poDS == NULL ){
    std::cerr<<"Failed to open '"<<filename<<"'!"<<std::endl;
    throw std::runtime_error("Failed to open shapefile!");
  }

  Polygons geometries;

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

  return geometries;
}

TEST_CASE("Shapefile open failure"){
  CHECK_THROWS(ReadShapefile("asdfjasfdlkjasdf", "ajdsfajsfdklajsfd"));
}



GeographicLib::Geodesic GreatCircleGenerator::geod = GeographicLib::Geodesic(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f());

GreatCircleGenerator::GreatCircleGenerator(
  const Point2D &a,
  const Point2D &b,
  const int num_pts0
){
  gline = geod.InverseLine(
    a.y*RAD_TO_DEG,
    a.x*RAD_TO_DEG,
    b.y*RAD_TO_DEG,
    b.x*RAD_TO_DEG
  );

  num_pts = num_pts0;

  da = gline.Arc() / num_pts;
}

Point2D GreatCircleGenerator::operator()(int i) const {
  Point2D temp;
  gline.ArcPosition(i * da, temp.y, temp.x);
  temp.toRadians();
  return temp;
}
  
int GreatCircleGenerator::size() const {
  return num_pts;
}
