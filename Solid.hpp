#ifndef _Solid_hpp_
#define _Solid_hpp_

#include "Point.hpp"
#include "Orientation.hpp"
#include <cmath>
#include <array>
#include <vector>
#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>

const double IEL = std::atan(0.5); //Icosahedron equatorial latitudes
const double IES = 36*M_PI/180;    //Icosahedron equatorial spacing
const double phi = (1+std::sqrt(5))/2;

class SolidXYZ;
class SolidXY{
 public:
  static const int verts = 12;
  std::array<Point2D,verts> v = {{
    { 0    ,  M_PI/2},
    { 0    , -M_PI/2},
    {-5*IES,     IEL},
    {-4*IES,    -IEL},
    {-3*IES,     IEL},
    {-2*IES,    -IEL},
    {-1*IES,     IEL},
    { 0*IES,    -IEL},
    { 1*IES,     IEL},
    { 2*IES,    -IEL},
    { 3*IES,     IEL},
    { 4*IES,    -IEL}
  }};

  SolidXY() = default;
  SolidXY(const Orientation &o);
  SolidXY& rotate(const Point2D &p, double rtheta);
  SolidXY& rotate(double rlat, double rlon, double rtheta);
  SolidXY& rotateTheta(const double rtheta);
  SolidXY& toMercator();
  SolidXY& toRadians();
  std::vector<int> neighbors() const;
  double neighborDistance() const;
  SolidXYZ toXYZ(const double radius) const;

  template <class Archive>
  void serialize( Archive & ar ){
    ar(v);
  }
};



class SolidXYZ {
 private:

 public:
  SolidXYZ() = default; //TODO: Make private?
  friend class SolidXY;

  //These values are a direct translation of those for SolidXY
  std::array<Point3D,SolidXY::verts> v = {{
    {std::nan(""), std::nan(""), std::nan("")},
    {std::nan(""), std::nan(""), std::nan("")},
    {std::nan(""), std::nan(""), std::nan("")},
    {std::nan(""), std::nan(""), std::nan("")},
    {std::nan(""), std::nan(""), std::nan("")},
    {std::nan(""), std::nan(""), std::nan("")},
    {std::nan(""), std::nan(""), std::nan("")},
    {std::nan(""), std::nan(""), std::nan("")},
    {std::nan(""), std::nan(""), std::nan("")},
    {std::nan(""), std::nan(""), std::nan("")},
    {std::nan(""), std::nan(""), std::nan("")},
    {std::nan(""), std::nan(""), std::nan("")}
  }};

  SolidXY toLatLon() const;
  SolidXYZ& rotateTo(const Point3D &o);
  std::vector<int> neighbors() const;

  template <class Archive>
  void serialize( Archive & ar ) {
    ar(v);
  }
};


#endif