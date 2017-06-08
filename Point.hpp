#ifndef _point_hpp_
#define _point_hpp_

#include <cereal/archives/binary.hpp>

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
  template <class Archive>
  void serialize( Archive & ar ) {
    ar(x,y);
  }
};

class Point3D {
 public:
  double x,y,z;
  Point3D(double x, double y, double z);
  Point2D toLatLon() const;
  template <class Archive>
  void serialize( Archive & ar ){
    ar(x,y,z);
  }
};

#endif