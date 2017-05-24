#ifndef _geostuff_hpp_
#define _geostuff_hpp_

#include "Point.hpp"

double GeoDistanceFlatEarth(const Point2D &a, const Point2D &b);

double GeoDistanceHaversine(const Point2D &pa, const Point2D &pb);

Point2D WGS84toEPSG3857(const Point2D &p);

template<class T>
void ToMercator(T &pts);

double EuclideanDistance(const Point3D &a, const Point3D &b);

#endif