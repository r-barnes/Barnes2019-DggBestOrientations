#ifndef _icosa_hpp_
#define _icosa_hpp_

#include "Point.hpp"
#include "POI.hpp"
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

  IcosaXY() = default;
  IcosaXY(const Orientation &o);
  IcosaXY& rotate(const Point2D &p, double rtheta);
  IcosaXY& rotate(double rlat, double rlon, double rtheta);
  IcosaXY& rotateTheta(const double rtheta);
  IcosaXY& toMercator();
  IcosaXY& toRadians();
  void print() const;
  std::vector<int> neighbors() const;
  double neighborDistance() const;
  IcosaXYZ toXYZ(const double radius) const;

  template <class Archive>
  void serialize( Archive & ar ){
    ar(v);
  }
};



class IcosaXYZ {
 private:

 public:
  IcosaXYZ() = default; //TODO: Make private?
  friend class IcosaXY;

  //These values are a direct translation of those for IcosaXY
  std::array<Point3D,IcosaXY::verts> v = {{
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

  IcosaXY toLatLon() const;
  void print() const;
  IcosaXYZ& rotateTo(const Point3D &o);
  std::vector<int> neighbors() const;

  template <class Archive>
  void serialize( Archive & ar ) {
    ar(v);
  }
};


#endif