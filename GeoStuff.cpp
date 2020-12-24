#include "GeoStuff.hpp"
#include <cmath>
#include "doctest.h"
#include <ogrsf_frmts.h>
#include <cmath>
#include <iostream>
#include <cassert>
#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/Constants.hpp>

const double RAD_TO_DEG = 180.0/M_PI;
const double Rearth = 6371; //km

// """
// Calculate the distance between two latitude-longitude points assuming a
// locally-approximated flat Earth.
// :param a  First point
// :param b  Second point
// :returns Distance between the two points in kilometres
// """
double GeoDistanceFlatEarth(
  const Point2D &a,
  const Point2D &b
){
  //Flat Earth Approx
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



// """
// Calculate the distance between two latitude-longitude points using the
// Haversine approximation.
// :param a  First point
// :param b  Second point
// :returns Distance between the two points in kilometres
// """
double GeoDistanceHaversine(
  const Point2D &pa,
  const Point2D &pb
){
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



// """
// Calculate the distance between two latitude-longitude points using a
// spherical approximation.
// :param a  First point
// :param b  Second point
// :returns Distance between the two points in kilometres
// """
double GeoDistanceSphere(
  const Point2D &a,
  const Point2D &b
){
  static auto geod = GeographicLib::Geodesic(1000*Rearth, 0);

  const auto gline = geod.InverseLine(
    a.y*RAD_TO_DEG,
    a.x*RAD_TO_DEG,
    b.y*RAD_TO_DEG,
    b.x*RAD_TO_DEG
  );

  return gline.Distance()/1000; //km
}



// """
// Calculate the distance between two latitude-longitude points using a
// high-precision ellipsoidal approximation.
// :param a  First point
// :param b  Second point
// :returns Distance between the two points in kilometres
// """
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
// """
// Convert a WGS84 lat-long point to a EPSG3857 Mercator x-y point
// :param p  Point to be converted
// :returns  A EPSG3857 Mercator x-y point
// """
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



// """
// Get the planar distance between two points in their units
// See: https://github.com/proj4js/proj4js/blob/2006b0a06d000308caa3625005f3d5734ef11f61/lib/projections/merc.js
// :param a  First point
// :param b  Second point
// :returns  Distance between the points
// """
double EuclideanDistance(const Point3D &a, const Point3D &b){
  return std::sqrt(std::pow(a.x-b.x,2)+std::pow(a.y-b.y,2)+std::pow(a.z-b.z,2));
}

TEST_CASE("EuclideanDistance"){
  auto a = Point3D(10, 33, 17);
  auto b = Point3D(97, 42, 3);
  auto dist = EuclideanDistance(a,b);
  CHECK(dist == doctest::Approx(88.5776495));
}



// """
// Convert a WGS84 lat-long point to an x-y-z Ellipsoidal Cartesian point
// :param p  Point to be converted
// :returns  An x-y-z point
// """
Point3D WGS84toEllipsoidCartesian(const Point2D &p) {
  static const GeographicLib::Geocentric& earth = GeographicLib::Geocentric::WGS84();

  Point3D temp;

  //earth.Forward(lat, lon, h, X, Y, Z);
  earth.Forward(p.y*RAD_TO_DEG,p.x*RAD_TO_DEG,0,temp.x,temp.y,temp.z);

  return temp;
}



// """
// Convert an x-y-z Ellipsoidal Cartesian point to a WGS84 lat-long point
// :param p  Point to be converted
// :returns  A WGS84 lat-long point
// """
Point2D EllipsoidCartesiantoWGS84(const Point3D &p){
  static const GeographicLib::Geocentric& earth = GeographicLib::Geocentric::WGS84();

  Point2D temp;
  double height;
  earth.Reverse(p.x,p.y,p.z,temp.y,temp.x,height);
  CHECK(height==doctest::Approx(0));
  temp.toRadians();
  return temp;
}



// """
// Convert a WGS84 lat-long point to a x-y-z spherial Cartesian point
// :param p  Point to be converted
// :returns  A WGS84 lat-long point
// """
Point3D WGS84toSphericalCartesian(const Point2D &p){
  return p.toXYZ(Rearth);
}



// """
// Convert an x-y-z spherical Cartesian point to a WGS84 lat-long point
// :param p  Point to be converted
// :returns  A WGS84 lat-long point
// """
Point2D SphericalCartesiantoWGS84(const Point3D &p){
  return p.toLatLon();
}



// """
// Reads in a shapefile and returns a set of polygons described by the file.
// This is useful for reading in coastline data.
// :param   filename  Name of file to be read
// :param   layername Layer to be read from the file
// :returns Vector of polygons from the file.
// """
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
    } else if( wkbFlatten(poGeometry->getGeometryType()) == wkbLineString ){
      OGRLineString *ls = (OGRLineString *) poGeometry;
      geometries.emplace_back();
      for(int i=0;i<ls->getNumPoints();i++)
        geometries.back().exterior.emplace_back(ls->getX(i),ls->getY(i));
    } else if( wkbFlatten(poGeometry->getGeometryType()) == wkbMultiLineString ){
      OGRMultiLineString *mls = (OGRMultiLineString *) poGeometry;
      for(const auto &ls: *mls){
        geometries.emplace_back();
        for(int i=0;i<ls->getNumPoints();i++)
          geometries.back().exterior.emplace_back(ls->getX(i),ls->getY(i));
      }
    } else if( wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon ){
      OGRPolygon *poly = (OGRPolygon *) poGeometry;
      auto extring = poly->getExteriorRing();
      //Ignore interior rings for now: they're probably lakes
      geometries.emplace_back();
      for(int i=0;i<extring->getNumPoints();i++)
        geometries.back().exterior.emplace_back(extring->getX(i),extring->getY(i));
    } else {
      std::cerr<<"Unrecognised geometry of type: "<<OGRGeometryTypeToName(poGeometry->getGeometryType())<<std::endl;
    }
    OGRFeature::DestroyFeature( poFeature );
  }
  GDALClose( poDS );

  return geometries;
}

TEST_CASE("Shapefile open failure"){
  CHECK_THROWS(ReadShapefile("asdfjasfdlkjasdf", "ajdsfajsfdklajsfd"));
}



// """
// :returns Number of points to be generated by the generator
// """
unsigned int GreatCircleGenerator::size() const {
  return num_pts;
}



// """
// :returns Spacing between the generator's points
// """
double GreatCircleGenerator::getSpacing() const {
  return spacing;
}



// """
// Construct a great circle generator
// :param  geod_radius     Radius of the geod on which the points are generated (see GeographicLib)
// :param  geod_flattening Flattening of the geod on which the points are generated (see GeographicLib)
// :param  a               Start point of the great circle
// :param  b               End point of the great circle
// :param  spacing0        Distance between interpolated points in metres
// """
GeographicLibGreatCircleGenerator::GeographicLibGreatCircleGenerator(
  double geod_radius,
  double geod_flattening,
  const Point2D &a,
  const Point2D &b,
  const double spacing0
){
  geod = GeographicLib::Geodesic(geod_radius, geod_flattening);

  gline = geod.InverseLine(
    a.y*RAD_TO_DEG,
    a.x*RAD_TO_DEG,
    b.y*RAD_TO_DEG,
    b.x*RAD_TO_DEG
  );

  spacing = spacing0;

  const auto dist = gline.Distance()/1000; //km
  assert(!std::isnan(dist));
  num_pts         = (int)std::ceil(dist/spacing);
  da              = gline.Arc()/num_pts;
}



// """
// Get ith the interpolated point on the great circle
// :param   i        Index of interpolated point to get
// :returns An interpolated point lying on the specified great circle
// """
Point2D GeographicLibGreatCircleGenerator::operator()(const int i) const {
  Point2D temp;
  gline.ArcPosition(i * da, temp.y, temp.x);
  temp.toRadians();
  return temp;
}



// """
// Construct an ellipsoidal great circle generator
// :param  a               Start point of the great circle
// :param  b               End point of the great circle
// :param  spacing0        Distance between interpolated points in metres
// """
EllipsoidalGreatCircleGenerator::EllipsoidalGreatCircleGenerator(const Point2D &a, const Point2D &b, const double spacing0) :
  GeographicLibGreatCircleGenerator(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f(), a, b, spacing0) {}



// """
// Construct a spherical great circle generator
// :param  a               Start point of the great circle
// :param  b               End point of the great circle
// :param  spacing0        Distance between interpolated points in metres
// """
SphericalGreatCircleGenerator::SphericalGreatCircleGenerator(const Point2D &a, const Point2D &b, const double spacing0) :
  GeographicLibGreatCircleGenerator(1000*Rearth, 0, a, b, spacing0) {}



// """
// Construct a spherical great circle generator without using GeographicLib
// :param  a               Start point of the great circle
// :param  b               End point of the great circle
// :param  spacing0        Distance between interpolated points in metres
// """
SimpleSphericalGreatCircleGenerator::SimpleSphericalGreatCircleGenerator(
  const Point2D &a,
  const Point2D &b,
  const double spacing0
){
  spacing = spacing0;

  svec = a.toXYZ(1);
  const Point3D n2 = b.toXYZ(1);

  const auto dist  = GeoDistanceHaversine(a,b); //TODO: Conflict this returns in kilometres
  num_pts          = std::floor(dist/spacing0)+1;

  const double sd  = svec.cross(n2).mag();
  const double cd  = svec.dot(n2);
  const double ang = std::atan2(sd,cd);
  da               = ang/num_pts;

  bvec = svec.cross(n2).unitify().cross(svec);

  CHECK(svec.mag2()!=doctest::Approx(0));
  CHECK(bvec.mag2()!=doctest::Approx(0));
};



// """
// Get the ith great circle interpolant without using GeographicLib
// :param  i        Index of the point to get
// :returns A point lying on the specified great circle
// """
Point2D SimpleSphericalGreatCircleGenerator::operator()(const int i) const {
  const double sd = std::sin(i*da);
  const double cd = std::cos(i*da);

  Point3D temp;
  temp.x = cd*svec.x + sd*bvec.x;
  temp.y = cd*svec.y + sd*bvec.y;
  temp.z = cd*svec.z + sd*bvec.z;

  CHECK(temp.mag2()!=doctest::Approx(0));

  temp.unitify();

  return temp.toLatLon();
}



// """
// Return a great circle generator of the specified type
// :param  gcgt        The type of great circle generator to make
// :param  a               Start point of the great circle
// :param  b               End point of the great circle
// :param  spacing0        Distance between interpolated points in metres
// :return A shared pointer to an appropriately initialized great circle generator
// """
std::shared_ptr<GreatCircleGenerator> GreatArcFactory::make(
  const GreatCircleGeneratorType gcgt,
  const Point2D &a,
  const Point2D &b,
  const double spacing0  //TODO: Is this in kilometres or meters? Because MAX_COAST_INTERPOINT_DIST is listed as being in kilometres
){
  if(gcgt == GreatCircleGeneratorType::SPHERICAL)
    return std::make_shared<SphericalGreatCircleGenerator>(a,b,spacing0);
  else if(gcgt == GreatCircleGeneratorType::ELLIPSOIDAL)
    return std::make_shared<EllipsoidalGreatCircleGenerator>(a,b,spacing0);
  else if(gcgt == GreatCircleGeneratorType::SIMPLE_SPHERICAL)
    return std::make_shared<SimpleSphericalGreatCircleGenerator>(a,b,spacing0);
  else
    throw std::runtime_error("Unrecognised projection!");
  return NULL;
}
