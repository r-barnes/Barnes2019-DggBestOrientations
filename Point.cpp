#include "Point.hpp"
#include <cmath>

const double DEG_TO_RAD = M_PI/180.0;
const double RAD_TO_DEG = 180.0/M_PI;

Point2D::Point2D(){}

Point2D::Point2D(double x, double y) {
  this->x = x;
  this->y = y;
}

void Point2D::toRadians() {
  x *= DEG_TO_RAD;
  y *= DEG_TO_RAD;
}

void Point2D::toDegrees() {
  x *= RAD_TO_DEG;
  y *= RAD_TO_DEG;
}

Point3D Point2D::toXYZ(const double radius) const {
  return Point3D(
    radius * std::cos(x) * std::cos(y),
    radius * std::sin(x) * std::cos(y),
    radius * std::sin(y)
  );
}

Point3D::Point3D(double x, double y, double z) {
  this->x = x;
  this->y = y;
  this->z = z;
}

Point2D Point3D::toLatLon() const {
  double radius = x*x + y*y + z*z; //TODO: Square root?
  return Point2D(
    std::atan2(y,x),
    std::asin(z/radius)
  );
}

