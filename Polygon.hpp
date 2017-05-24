#ifndef _Polygon_hpp_
#define _Polygon_hpp_

#include <vector>
#include <limits>
#include "Point.hpp"

class Polygon {
 public:
  std::vector<Point2D> exterior;
  double minX() const;
  double maxX() const;
  double minY() const;
  double maxY() const;
  void toRadians();
  void toDegrees();
  bool containsPoint(const Point2D &xy) const;
};

typedef std::vector<Polygon> Polygons;

#endif