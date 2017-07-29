#ifndef _point_hpp_
#define _point_hpp_

class Point3D;

class Point2D {
 public:
  double x;
  double y;
  Point2D() = default;
  Point2D(double x0, double y0);
  Point2D& toRadians();
  Point2D& toDegrees();
  Point2D& rotateTheta(const double rtheta);
  Point3D toXYZ(const double radius) const;
};

class Point3D {
 public:
  double x,y,z;
  Point3D() = default;
  Point3D(double x0, double y0, double z0);
  Point2D toLatLon() const;
  double dot(const Point3D &b) const;
};

class Rotator {
 private:
  double mr;
  double R_a;
  double R_b;
  double R_c;
  double R_d;
  double R_e;
  double R_f;
  double R_g;
  double R_h;
  double R_i;
 public:
  Rotator(const Point3D &oldv, const Point3D &newv);
  Point3D operator()(const Point3D &p) const;
};

#endif