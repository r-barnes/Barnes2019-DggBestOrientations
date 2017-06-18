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

IcosaXY::IcosaXY(const Orientation &o){
  rotate(o.pole,o.theta);
}



IcosaXY& IcosaXY::rotate(const Point2D &p, double rtheta){
  rotateTheta(rtheta);
  *this = toXYZ(1).rotateTo(p.toXYZ(1)).toLatLon();
  return *this;
}

TEST_CASE("rotate"){
  const auto a = IcosaXY();
  const auto p = IcosaXY().rotate(Point2D(0,90).toRadians(),0);
  for(unsigned int i=0;i<a.v.size();i++){
    CHECK(a.v[i].x==doctest::Approx(p.v[i].x));
    CHECK(a.v[i].y==doctest::Approx(p.v[i].y));
  }
}



IcosaXY& IcosaXY::rotate(double rlat, double rlon, double rtheta){
  Point2D temp(rlon,rlat);
  return rotate(temp,rtheta);
}



IcosaXY& IcosaXY::rotateTheta(const double rtheta){
  if(rtheta==0)
    return *this;
  
  assert(std::abs(v[0].y-M_PI/2)<1e-6); //Can only rotate when icosahedron is North-South aligned

  for(auto &p: v)
    p.rotateTheta(rtheta);
  return *this;
}

TEST_CASE("rotateTheta"){
  IcosaXY().rotateTheta(23*DEG_TO_RAD);
}



IcosaXY& IcosaXY::toMercator(){
  for(auto &p: v)
    p = WGS84toEPSG3857(p);
  return *this;
}

TEST_CASE("toMercator"){
  auto ico2d = IcosaXY();
  auto p     = WGS84toEPSG3857(ico2d.v[3]);
  ico2d.toMercator();
  CHECK(p.x==ico2d.v[3].x);
  CHECK(p.y==ico2d.v[3].y);
}



IcosaXY& IcosaXY::toRadians(){
  for(auto &p: v)
    p.toRadians();
  return *this;
}

TEST_CASE("toRadians"){
  auto ico2d = IcosaXY();
  auto p     = ico2d.v[3];
  ico2d.toRadians();
  p.toRadians();
  CHECK(p.x==ico2d.v[3].x);
  CHECK(p.y==ico2d.v[3].y);
}



void IcosaXY::print() const {
  for(const auto &p: v)
    std::cerr<<std::fixed<<std::setw(10)<<p.y<<" "<<std::fixed<<std::setw(10)<<p.x<<" -- "<<std::setw(10)<<std::fixed<<p.y*RAD_TO_DEG<<" "<<std::setw(10)<<std::fixed<<p.x*RAD_TO_DEG<<std::endl;
}



//Return neighbours in a N1a,N1b,N2a,N2b,... format
std::vector<int> IcosaXY::neighbors() const {
  std::vector<int> ret;
  double dist = std::numeric_limits<double>::infinity();
  //Find nearest vertex that isn't itself
  for(unsigned int i=1;i<v.size();i++)
    dist = std::min(dist,GeoDistanceHaversine(v[0],v[i]));
  //Increase by 5% to account for floating point errors
  dist *= 1.05;
  for(unsigned int i=0;  i<v.size();i++)
  for(unsigned int j=i+1;j<v.size();j++)
    if(GeoDistanceHaversine(v[i],v[j])<dist){ //IcosaXY is about as close as any close pole
      ret.emplace_back(i);
      ret.emplace_back(j);
    }
  return ret;
}

TEST_CASE("neighbors"){
  const auto n = IcosaXY().neighbors();
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



double IcosaXY::neighborDistance() const {
  const auto n = neighbors();
  return GeoDistanceHaversine(v[n[0]], v[n[1]]);
}

TEST_CASE("neighborDistance"){
  IcosaXY icoxy;
  const auto n     = icoxy.neighbors();
  const auto ndist = icoxy.neighborDistance();
  for(unsigned int i=0;i<n.size();i+=2)
    CHECK(GeoDistanceHaversine(icoxy.v[n.at(i)],icoxy.v[n.at(i+1)])==doctest::Approx(ndist));
}



IcosaXYZ IcosaXY::toXYZ(const double radius) const {
  IcosaXYZ temp;
  for(unsigned int i=0;i<v.size();i++)
    temp.v[i] = v[i].toXYZ(radius);
  return temp;
}

TEST_CASE("toXYZ"){
  auto a = IcosaXY();
  auto p = IcosaXY().toXYZ(1).toLatLon();
  for(unsigned int i=0;i<a.v.size();i++){
    CHECK(a.v[i].x==doctest::Approx(p.v[i].x));
    CHECK(a.v[i].y==doctest::Approx(p.v[i].y));
  }
}



IcosaXY IcosaXYZ::toLatLon() const {
  IcosaXY temp;
  for(unsigned int i=0;i<v.size();i++)
    temp.v[i] = v[i].toLatLon();
  return temp;
}

TEST_CASE("IcosaXYZ::toLatLon"){
  const auto a = IcosaXY();
  const auto p = IcosaXY().toXYZ(1).toLatLon();
  for(unsigned int i=0;i<a.v.size();i++){
    CHECK(a.v[i].x==doctest::Approx(p.v[i].x));
    CHECK(a.v[i].y==doctest::Approx(p.v[i].y));
  }
}



void IcosaXYZ::print() const {
  for(const auto &p: v)
    std::cerr<<std::fixed<<std::setprecision(10)<<std::setw(15)<<p.x<<" "
             <<std::fixed<<std::setprecision(10)<<std::setw(15)<<p.y<<" "
             <<std::fixed<<std::setprecision(10)<<std::setw(15)<<p.z
             <<std::endl;
}

//https://math.stackexchange.com/a/476311/14493
IcosaXYZ& IcosaXYZ::rotateTo(const Point3D &o){
  Rotator r(v[0],o);

  //Rotate each pole
  for(auto &p: v)
    p = r(p);

  return *this;
}

TEST_CASE("rotateTo"){
  SUBCASE("theta rotation"){
    const auto a = IcosaXY();
    auto r = IcosaXY().rotate(90*DEG_TO_RAD,0,45); //Rotate forward
    r.rotate(90*DEG_TO_RAD,0,-45);                 //Rotate back
    for(unsigned int i=0;i<a.v.size();i++){
      CHECK(a.v[i].x==doctest::Approx(r.v[i].x));
      CHECK(a.v[i].y==doctest::Approx(r.v[i].y));
    }
  }
}



std::vector<int> IcosaXYZ::neighbors() const {
  std::vector<int> ret;
  double dist = std::numeric_limits<double>::infinity();
  //Find nearest vertex that isn't itself
  for(unsigned int i=1;i<v.size();i++)
    dist = std::min(dist,EuclideanDistance(v[0],v[i]));
  //Increase by 10% to account for floating point errors
  dist *= 1.1;
  for(unsigned int i=0;  i<v.size();i++)
  for(unsigned int j=i+1;j<v.size();j++)
    if(EuclideanDistance(v[i],v[j])<dist){ //IcosaXY is about as close as any close pole
      ret.emplace_back(i);
      ret.emplace_back(j);
    }
  return ret;
}



TEST_CASE("neighbors"){
  const auto n = IcosaXY().toXYZ(1).neighbors();
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
