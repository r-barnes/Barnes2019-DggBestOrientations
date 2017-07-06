#ifndef _geostuff_hpp_
#define _geostuff_hpp_

#include "Point.hpp"
#include "Polygon.hpp"
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicLine.hpp>

double GeoDistanceFlatEarth(const Point2D &a, const Point2D &b);

double GeoDistanceHaversine(const Point2D &pa, const Point2D &pb);

double GeoDistanceEllipsoid(const Point2D &a, const Point2D &b);

Point2D WGS84toEPSG3857(const Point2D &p);

Point3D WGS84toEllipsoidCartesian(const Point2D &p);

Point2D EllipsoidCartesiantoWGS84(const Point3D &p);

template<class T>
void ToMercator(T &pts);

double EuclideanDistance(const Point3D &a, const Point3D &b);

Polygons ReadShapefile(std::string filename, std::string layername);

class GreatCircleGenerator {
 public:
  static GeographicLib::Geodesic geod;
 private:
  GeographicLib::GeodesicLine gline;
  double da;
  int num_pts;
 public:
  GreatCircleGenerator(const Point2D &a, const Point2D &b, const int num_pts0);
  Point2D operator()(int i) const;
  int size() const;
};

#endif