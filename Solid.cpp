#include "Solid.hpp"
#include "GeoStuff.hpp"
#include <iomanip>
#include <iostream>
#include <limits>
#include <cassert>
#include "doctest.h"

#include <cereal/types/array.hpp>

const double DEG_TO_RAD = M_PI/180.0;
const double RAD_TO_DEG = 180.0/M_PI;

SolidXY::SolidXY(const Orientation &o){
  rotate(o.pole,o.theta);
}



SolidXY& SolidXY::rotate(const Point2D &p, double rtheta){
  rotateTheta(rtheta);
  *this = toXYZ(1).rotateTo(p.toXYZ(1)).toLatLon();
  return *this;
}

TEST_CASE("rotate"){
  const auto a = SolidXY();
  const auto p = SolidXY().rotate(Point2D(0,90).toRadians(),0);
  for(unsigned int i=0;i<a.v.size();i++){
    CHECK(a.v[i].x==doctest::Approx(p.v[i].x));
    CHECK(a.v[i].y==doctest::Approx(p.v[i].y));
  }
}



SolidXY& SolidXY::rotate(double rlat, double rlon, double rtheta){
  Point2D temp(rlon,rlat);
  return rotate(temp,rtheta);
}



SolidXY& SolidXY::rotateTheta(const double rtheta){
  if(rtheta==0)
    return *this;
  
  assert(std::abs(v[0].y-M_PI/2)<1e-6); //Can only rotate when icosahedron is North-South aligned

  for(auto &p: v)
    p.rotateTheta(rtheta);
  return *this;
}

TEST_CASE("rotateTheta"){
  SolidXY().rotateTheta(23*DEG_TO_RAD);
}



SolidXY& SolidXY::toMercator(){
  for(auto &p: v)
    p = WGS84toEPSG3857(p);
  return *this;
}

TEST_CASE("toMercator"){
  auto ico2d = SolidXY();
  auto p     = WGS84toEPSG3857(ico2d.v[3]);
  ico2d.toMercator();
  CHECK(p.x==ico2d.v[3].x);
  CHECK(p.y==ico2d.v[3].y);
}



SolidXY& SolidXY::toRadians(){
  for(auto &p: v)
    p.toRadians();
  return *this;
}

TEST_CASE("toRadians"){
  auto ico2d = SolidXY();
  auto p     = ico2d.v[3];
  ico2d.toRadians();
  p.toRadians();
  CHECK(p.x==ico2d.v[3].x);
  CHECK(p.y==ico2d.v[3].y);
}



//Return neighbours in a N1a,N1b,N2a,N2b,... format
std::vector<int> SolidXY::neighbors() const {
  std::vector<int> ret;
  double dist = std::numeric_limits<double>::infinity();
  //Find nearest vertex that isn't itself
  for(unsigned int i=1;i<v.size();i++)
    dist = std::min(dist,GeoDistanceHaversine(v[0],v[i]));
  //Increase by 5% to account for floating point errors
  dist *= 1.05;
  for(unsigned int i=0;  i<v.size();i++)
  for(unsigned int j=i+1;j<v.size();j++)
    if(GeoDistanceHaversine(v[i],v[j])<dist){ //SolidXY is about as close as any close pole
      ret.emplace_back(i);
      ret.emplace_back(j);
    }
  return ret;
}

TEST_CASE("neighbors"){
  const auto n = SolidXY().neighbors();
  CHECK(n.size()==2*30); //I should have 30 edges represented by 60 neighbour pairs
  
  std::vector< std::vector<int> > ncheck(12);
  for(unsigned int ni=0;ni<n.size();ni+=2){
    ncheck.at(n.at(ni)).emplace_back(n.at(ni+1));
    ncheck.at(n.at(ni+1)).emplace_back(n.at(ni));
  }

  //Each vertex should have 5 neighbours
  for(const auto &ni: ncheck)
    CHECK(ni.size()==5);
}



double SolidXY::neighborDistance() const {
  const auto n = neighbors();
  return GeoDistanceHaversine(v[n[0]], v[n[1]]);
}

TEST_CASE("neighborDistance"){
  SolidXY icoxy;
  const auto n     = icoxy.neighbors();
  const auto ndist = icoxy.neighborDistance();
  for(unsigned int i=0;i<n.size();i+=2)
    CHECK(GeoDistanceHaversine(icoxy.v[n.at(i)],icoxy.v[n.at(i+1)])==doctest::Approx(ndist));
}



SolidXYZ SolidXY::toXYZ(const double radius) const {
  SolidXYZ temp;
  for(unsigned int i=0;i<v.size();i++)
    temp.v[i] = v[i].toXYZ(radius);
  return temp;
}

TEST_CASE("toXYZ"){
  auto a = SolidXY();
  auto p = SolidXY().toXYZ(1).toLatLon();
  for(unsigned int i=0;i<a.v.size();i++){
    CHECK(a.v[i].x==doctest::Approx(p.v[i].x));
    CHECK(a.v[i].y==doctest::Approx(p.v[i].y));
  }
}



SolidXY SolidXYZ::toLatLon() const {
  SolidXY temp;
  for(unsigned int i=0;i<v.size();i++)
    temp.v[i] = v[i].toLatLon();
  return temp;
}

TEST_CASE("SolidXYZ::toLatLon"){
  const auto a = SolidXY();
  const auto p = SolidXY().toXYZ(1).toLatLon();
  for(unsigned int i=0;i<a.v.size();i++){
    CHECK(a.v[i].x==doctest::Approx(p.v[i].x));
    CHECK(a.v[i].y==doctest::Approx(p.v[i].y));
  }
}



//https://math.stackexchange.com/a/476311/14493
SolidXYZ& SolidXYZ::rotateTo(const Point3D &o){
  Rotator r(v[0],o);

  //Rotate each pole
  for(auto &p: v)
    p = r(p);

  return *this;
}

TEST_CASE("rotateTo"){
  SUBCASE("theta rotation"){
    const auto a = SolidXY();
    auto r = SolidXY().rotate(90*DEG_TO_RAD,0,45); //Rotate forward
    r.rotate(90*DEG_TO_RAD,0,-45);                 //Rotate back
    for(unsigned int i=0;i<a.v.size();i++){
      CHECK(a.v[i].x==doctest::Approx(r.v[i].x));
      CHECK(a.v[i].y==doctest::Approx(r.v[i].y));
    }
  }
}



std::vector<int> SolidXYZ::neighbors() const {
  std::vector<int> ret;
  double dist = std::numeric_limits<double>::infinity();
  //Find nearest vertex that isn't itself
  for(unsigned int i=1;i<v.size();i++)
    dist = std::min(dist,EuclideanDistance(v[0],v[i]));
  //Increase by 10% to account for floating point errors
  dist *= 1.1;
  for(unsigned int i=0;  i<v.size();i++)
  for(unsigned int j=i+1;j<v.size();j++)
    if(EuclideanDistance(v[i],v[j])<dist){ //SolidXY is about as close as any close pole
      ret.emplace_back(i);
      ret.emplace_back(j);
    }
  return ret;
}



TEST_CASE("neighbors"){
  const auto n = SolidXY().toXYZ(1).neighbors();
  CHECK(n.size()==2*30); //I should have 30 edges represented by 60 neighbour pairs
  
  std::vector< std::vector<int> > ncheck(12);
  for(unsigned int ni=0;ni<n.size();ni+=2){
    ncheck.at(n.at(ni)).emplace_back(n.at(ni+1));
    ncheck.at(n.at(ni+1)).emplace_back(n.at(ni));
  }

  //Each vertex should have 5 neighbours
  for(const auto &ni: ncheck)
    CHECK(ni.size()==5);
}
