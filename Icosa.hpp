#ifndef _icosa_hpp_
#define _icosa_hpp_

#include "Point.hpp"
#include <cmath>
#include <array>
#include <vector>
#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>

const double IEL = std::atan(0.5); //Icosahedron equatorial latitudes
const double IES = 36*M_PI/180;    //Icosahedron equatorial spacing
const double phi = (1+std::sqrt(5))/2;

class IcosaXYZ;
class IcosaXY{
 public:
  std::array<Point2D,12> v = {{
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

  IcosaXY();
  IcosaXY(double rlat, double rlon, double rtheta);
  IcosaXY(const Point2D &p, double rtheta);
  IcosaXY& rotate(const Point2D &p, double rtheta);
  IcosaXY& rotate(double rlat, double rlon, double rtheta);
  IcosaXY& rotateTheta(const double rtheta);
  void toMercator();
  void toRadians();
  void print() const;
  std::vector<int> neighbors() const;
  double neighborDistance() const;
  IcosaXYZ toXYZ() const;

  template <class Archive>
  void serialize( Archive & ar ){
    ar(v);
  }
};

class IcosaXYZ {
 public:
  //These values are a direct translation of those for IcosaXY
  std::array<Point3D,12> v = {{
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

  IcosaXYZ();
  IcosaXYZ(double rlat, double rlon, double rtheta);
  IcosaXY toLatLon() const;
  void print() const;
  IcosaXYZ& rotateTo(const Point3D &o);
  void rotate(double rlat, double rlon, double rtheta);
  std::vector<int> neighbors() const;

  template <class Archive>
  void serialize( Archive & ar ) {
    ar(v);
  }
};


#endif