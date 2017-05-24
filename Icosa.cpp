#include "Icosa.hpp"
#include "GeoStuff.hpp"
#include <iomanip>
#include <iostream>
#include <limits>

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

void IcosaXY::rotate(double rlat, double rlon, double rtheta){
  for(auto &pole: v)
    pole = RotatePoint(rlat, rlon, rtheta, pole);
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
  //Increase by 10% to account for floating point errors
  dist *= 1.1;
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

//https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
IcosaXYZ& IcosaXYZ::rotateTo(const Point3D &o){
  const double x2 = v[0].x;
  const double y2 = v[0].y;
  const double z2 = v[0].z;
  const double x1 = o.x;
  const double y1 = o.y;
  const double z1 = o.z;
  const double a  = ((-x1*y2 + x2*y1)*(x1*y2 - x2*y1) + (-x1*z2 + x2*z1)*(x1*z2 - x2*z1))/(x1*x2 + y1*y2 + z1*z2 + 1) + 1;
  const double b  = -x1*y2 + x2*y1 + (-x1*z2 + x2*z1)*(y1*z2 - y2*z1)/(x1*x2 + y1*y2 + z1*z2 + 1);
  const double c  = -x1*z2 + x2*z1 + (-x1*y2 + x2*y1)*(-y1*z2 + y2*z1)/(x1*x2 + y1*y2 + z1*z2 + 1);
  const double d  = x1*y2 - x2*y1 + (x1*z2 - x2*z1)*(-y1*z2 + y2*z1)/(x1*x2 + y1*y2 + z1*z2 + 1);
  const double e  = ((-x1*y2 + x2*y1)*(x1*y2 - x2*y1) + (-y1*z2 + y2*z1)*(y1*z2 - y2*z1))/(x1*x2 + y1*y2 + z1*z2 + 1) + 1;
  const double f  = -y1*z2 + y2*z1 + (x1*y2 - x2*y1)*(-x1*z2 + x2*z1)/(x1*x2 + y1*y2 + z1*z2 + 1);
  const double g  = x1*z2 - x2*z1 + (x1*y2 - x2*y1)*(y1*z2 - y2*z1)/(x1*x2 + y1*y2 + z1*z2 + 1);
  const double h  = y1*z2 - y2*z1 + (-x1*y2 + x2*y1)*(x1*z2 - x2*z1)/(x1*x2 + y1*y2 + z1*z2 + 1);
  const double i  = ((-x1*z2 + x2*z1)*(x1*z2 - x2*z1) + (-y1*z2 + y2*z1)*(y1*z2 - y2*z1))/(x1*x2 + y1*y2 + z1*z2 + 1) + 1;
  for(auto &p: v)
    p = Point3D(
      a*p.x+b*p.y+c*p.z,
      d*p.x+e*p.y+f*p.z,
      g*p.x+h*p.y+i*p.z
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