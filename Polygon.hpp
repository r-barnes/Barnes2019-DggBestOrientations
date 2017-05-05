#ifndef _Polygon_hpp_
#define _Polygon_hpp_

#include <vector>
#include <limits>

class Point2D {
 public:
  double x;
  double y;
  Point2D(double x, double y);
};

class Polygon {
 public:
  std::vector<Point2D> exterior;
  double minX() const;
  double maxX() const;
  double minY() const;
  double maxY() const;
  void toRadians();
  void toDegrees();
  bool containsPoint(const double testx, const double testy) const;
  double distanceFromPoint(const double px, const double py) const;
};

#endif