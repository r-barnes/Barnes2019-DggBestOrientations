#include "Polygon.hpp"
#include <cmath>

const double DEG_TO_RAD = M_PI/180.0;
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

Point2D::Point2D(double x, double y) {
  this->x = x;
  this->y = y;
}

double Polygon::minX() const {
  double minx=std::numeric_limits<double>::infinity();
  for(const auto &p: exterior)
    minx = std::min(p.x,minx);
  return minx;
}

double Polygon::maxX() const {
  double maxx=-std::numeric_limits<double>::infinity();
  for(const auto &p: exterior)
    maxx = std::max(p.x,maxx);
  return maxx;
}

double Polygon::minY() const {
  double miny=std::numeric_limits<double>::infinity();
  for(const auto &p: exterior)
    miny=std::min(p.y,miny);
  return miny;
}

double Polygon::maxY() const {
  double maxy=-std::numeric_limits<double>::infinity();
  for(const auto &p: exterior)
    maxy=std::max(p.y,maxy);
  return maxy;
}

void Polygon::toRadians() {
  for(auto &p: exterior){
    p.x *= DEG_TO_RAD;
    p.y *= DEG_TO_RAD;
  }
}

void Polygon::toDegrees() {
  for(auto &p: exterior){
    p.x *= RAD_TO_DEG;
    p.y *= RAD_TO_DEG;
  }
}

//Info: https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
//Info: https://stackoverflow.com/questions/8721406/how-to-determine-if-a-point-is-inside-a-2d-convex-polygon/23223947#23223947
//Info: http://stackoverflow.com/a/2922778/752843
bool Polygon::containsPoint(const double testx, const double testy) const {
  unsigned int i, j;
  int c = 0;
  for (i = 0, j = exterior.size()-1; i < exterior.size(); j = i++) {
    if ( ((exterior[i].y>testy) != (exterior[j].y>testy)) &&
     (testx < (exterior[j].x-exterior[i].x) * (testy-exterior[i].y) / (exterior[j].y-exterior[i].y) + exterior[i].x) )
       c = !c;
  }
  return c;
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
double Polygon::distanceFromPoint(const double px, const double py) const {
  double dist = std::numeric_limits<double>::infinity();
  for(const auto &e: exterior)
    dist = std::min(dist,GeoDistance(px,py,e.x,e.y));
  return dist;
}
