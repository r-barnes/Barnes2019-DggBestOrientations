#ifndef _Solid_hpp_
#define _Solid_hpp_

#include "Point.hpp"
#include "Orientation.hpp"
#include <cmath>
#include <vector>
#include <functional>

class SolidXY;

SolidXY OrientationFullerIcosahedron();
SolidXY OrientationToRegularIcosahedron(const Orientation &o);
SolidXY OrientationToRegularDodecahedron(const Orientation &o);
SolidXY OrientationToRegularTetrahedron(const Orientation &o);
SolidXY OrientationToRegularOctahedron(const Orientation &o);
SolidXY OrientationToCuboctahedron(const Orientation &o);
SolidXY OrientationToPoint(const Orientation &o);

class SolidXYZ;
class SolidXY{
 private:
  friend class SolidXYZ;
  friend SolidXY OrientationFullerIcosahedron();
  friend SolidXY OrientationToPoint(const Orientation &o);
  friend SolidXY OrientationToRegularIcosahedron(const Orientation &o);
  friend SolidXY OrientationToRegularDodecahedron(const Orientation &o);
  friend SolidXY OrientationToRegularTetrahedron(const Orientation &o);
  friend SolidXY OrientationToRegularOctahedron(const Orientation &o);
  friend SolidXY OrientationToCuboctahedron(const Orientation &o);
  SolidXY() = default;

 public:
  std::vector<Point2D> v;

  SolidXY& rotate(const Point2D &p, double rtheta);
  SolidXY& rotate(double rlat, double rlon, double rtheta);
  SolidXY& rotateTheta(const double rtheta);
  SolidXY& toMercator();
  SolidXY& toRadians();
  std::vector<int> neighbors() const;
  double neighborDistance() const;
  SolidXYZ toXYZ(const double radius) const;
};



class SolidXYZ {
 private:

 public:
  SolidXYZ() = default; //TODO: Make private?
  friend class SolidXY;

  //These values are a direct translation of those for SolidXY
  std::vector<Point3D> v;

  SolidXY toLatLon() const;
  SolidXYZ& rotateTo(const Point3D &o);
  SolidXYZ& rotate(const Rotator &r);
  SolidXYZ& normalize();
  std::vector<int> neighbors() const;
};



typedef std::function<SolidXY(const Orientation &o)> SolidifyingFunc;

#endif