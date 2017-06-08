#include "Icosa.hpp"
#include "GeoStuff.hpp"
#include <iomanip>
#include <iostream>
#include <limits>

#include <cereal/types/array.hpp>

const double DEG_TO_RAD = M_PI/180.0;
const double RAD_TO_DEG = 180.0/M_PI;

//https://gis.stackexchange.com/questions/10808/lon-lat-transformation
Point2D RotatePoint(const double rlat, const double rlon, const double rtheta, Point2D ll){
  ll.x += M_PI+rtheta;                                       //Move [0,360] and add theta
  ll.x  = std::fmod(2*M_PI+std::fmod(ll.x,2*M_PI),2*M_PI); //Map back to [0,360]
  ll.x -= M_PI;                                              //Move back to [-180,180] system

  auto pr = ll.toXYZ(1);
  Point3D p(
     std::cos(rlat)*std::cos(rlon)*pr.x + std::sin(rlon)*pr.y + std::sin(rlat)*std::cos(rlon)*pr.z,
    -std::cos(rlat)*std::sin(rlon)*pr.x + std::cos(rlon)*pr.y - std::sin(rlat)*std::sin(rlon)*pr.z,
    -std::sin(rlat)*pr.x + std::cos(rlat)*pr.z
  );
  return p.toLatLon();
}

Point3D RotatePoint(const double rlat, const double rlon, const double rtheta, Point3D p3){
  return Point3D(
     std::cos(rlat)*std::cos(rlon)*p3.x + std::sin(rlon)*p3.y + std::sin(rlat)*std::cos(rlon)*p3.z,
    -std::cos(rlat)*std::sin(rlon)*p3.x + std::cos(rlon)*p3.y - std::sin(rlat)*std::sin(rlon)*p3.z,
    -std::sin(rlat)*p3.x + std::cos(rlat)*p3.z
  );
}

IcosaXY::IcosaXY(){}

IcosaXY::IcosaXY(double rlat, double rlon, double rtheta){
  rotate(rlat,rlon,rtheta);
}

IcosaXY::IcosaXY(const Point2D &p, double rtheta){
  rotate(p,rtheta);
}

IcosaXY& IcosaXY::rotate(const Point2D &p, double rtheta){
  rotateTheta(rtheta);
  *this = toXYZ().rotateTo(p.toXYZ(1)).toLatLon();
  return *this;
}

IcosaXY& IcosaXY::rotate(double rlat, double rlon, double rtheta){
  Point2D temp(rlon,rlat);
  return rotate(temp,rtheta);
}

IcosaXY& IcosaXY::rotateTheta(const double rtheta){
  for(auto &p: v){
    p.x += M_PI+rtheta;                                    //Move [0,360] and add theta
    p.x  = std::fmod(2*M_PI+std::fmod(p.x,2*M_PI),2*M_PI); //Map back to [0,360]
    p.x -= M_PI;                                           //Move back to [-180,180] system
  }
  return *this;
}

void IcosaXY::toMercator(){
  for(auto &p: v)
    p = WGS84toEPSG3857(p);
}

void IcosaXY::toRadians(){
  for(auto &p: v)
    p.toRadians();
}

void IcosaXY::print() const {
  for(const auto &p: v)
    std::cerr<<std::fixed<<std::setw(10)<<p.y<<" "<<std::fixed<<std::setw(10)<<p.x<<" -- "<<std::setw(10)<<std::fixed<<p.y*RAD_TO_DEG<<" "<<std::setw(10)<<std::fixed<<p.x*RAD_TO_DEG<<std::endl;
}

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

double IcosaXY::neighborDistance() const {
  auto n = neighbors();
  return GeoDistanceHaversine(v[n[0]], v[n[1]]);
}

IcosaXYZ IcosaXY::toXYZ() const {
  IcosaXYZ temp;
  for(unsigned int i=0;i<v.size();i++)
    temp.v[i] = v[i].toXYZ(1);
  return temp;
}





IcosaXYZ::IcosaXYZ(){}

IcosaXYZ::IcosaXYZ(double rlat, double rlon, double rtheta){
  rotate(rlat,rlon,rtheta);
}

IcosaXY IcosaXYZ::toLatLon() const {
  IcosaXY temp;
  for(unsigned int i=0;i<v.size();i++)
    temp.v[i] = v[i].toLatLon();
  return temp;
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
  //Icosahedron's North pole
  const double a1 = v[0].x;
  const double a2 = v[0].y;
  const double a3 = v[0].z;

  //Vector defining the new pole
  const double b1 = o.x;
  const double b2 = o.y;
  const double b3 = o.z;

  //B x A
  double r1 = -a2*b3 + a3*b2;
  double r2 =  a1*b3 - a3*b1;
  double r3 = -a1*b2 + a2*b1;
  const double mr = std::sqrt(r1*r1+r2*r2+r3*r3);

  if(mr==0)
    return *this;

  r1 /= mr;
  r2 /= mr;
  r3 /= mr;

  const double c = a1*b1 + a2*b2 + a3*b3;  //cos theta = B dot A
  const double s = sqrt(1-c*c);            //= sin theta
  const double t = 1-c;

  //Rotation Matrix from
  //Glassner, A.S. (Ed.), 1993. Graphics Gems I, 1st ed. p. 466
  const double R_a = c + r1*r1*t;
  const double R_b = r1*r2*t + r3*s;
  const double R_c = r1*r3*t - r2*s;
  const double R_d = r1*r2*t - r3*s;
  const double R_e = c + r2*r2*t;
  const double R_f = r1*s + r2*r3*t;
  const double R_g = r1*r3*t + r2*s;
  const double R_h = -r1*s + r2*r3*t;
  const double R_i = c + r3*r3*t;

  //Rotate each pole
  for(auto &p: v)
    p = Point3D(
      R_a*p.x+R_b*p.y+R_c*p.z,
      R_d*p.x+R_e*p.y+R_f*p.z,
      R_g*p.x+R_h*p.y+R_i*p.z
    );

  return *this;
}

void IcosaXYZ::rotate(double rlat, double rlon, double rtheta){
  for(auto &pole: v)
    pole = RotatePoint(rlat, rlon, rtheta, pole);
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