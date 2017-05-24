#ifndef _point_hpp_
#define _point_hpp_

class Point3D;

class Point2D {
 public:
  double x;
  double y;
  Point2D();
  Point2D(double x, double y);
  void toRadians();
  void toDegrees();
  Point3D toXYZ(const double radius) const;
};

class Point3D {
 public:
  double x,y,z;
  Point3D(double x, double y, double z);
  Point2D toLatLon() const;
};

#endif