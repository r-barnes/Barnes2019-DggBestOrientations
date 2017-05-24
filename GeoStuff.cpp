#include "GeoStuff.hpp"
#include <cmath>

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

double EuclideanDistance(const Point3D &a, const Point3D &b){
  return std::sqrt(std::pow(a.x-b.x,2)+std::pow(a.y-b.y,2)+std::pow(a.z-b.z,2));
}
