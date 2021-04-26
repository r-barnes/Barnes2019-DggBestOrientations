#include "Point.hpp"
#include <cmath>
#include <cassert>
#include <cstdlib>
#include "doctest.h"
#include <vector>
#include <stdexcept>

const double DEG_TO_RAD = M_PI/180.0;
const double RAD_TO_DEG = 180.0/M_PI;

Rotator::Rotator(const Point3D &oldv, const Point3D &newv){
  //Ensure that incoming vectors are normalized
  CHECK(oldv.mag2()==doctest::Approx(1));
  CHECK(newv.mag2()==doctest::Approx(1));

  Point3D r = oldv.cross(newv);

  mr = r.mag();

  if(mr==0)
    return;

  r.unitify();

  //Ensure that things were actually normalized
  CHECK(r.mag2()==doctest::Approx(1));

  const double c = oldv.dot(newv);   //cos theta = B dot A
  const double s = std::sqrt(1-c*c); //= sin theta
  const double t = 1-c;

  CHECK(s*s+c*c==doctest::Approx(1));

  //Rotation Matrix from
  //Glassner, A.S. (Ed.), 1993. Graphics Gems I, 1st ed. p. 466
  R_a = c + r.x*r.x*t;
  R_b = r.x*r.y*t + r.z*s;
  R_c = r.x*r.z*t - r.y*s;
  R_d = r.x*r.y*t - r.z*s;
  R_e = c + r.y*r.y*t;
  R_f = r.x*s + r.y*r.z*t;
  R_g = r.x*r.z*t + r.y*s;
  R_h = -r.x*s + r.y*r.z*t;
  R_i = c + r.z*r.z*t;
}

Point3D Rotator::operator()(const Point3D &p) const {
  if(mr==0)
    return p;

  CHECK(p.mag2()==doctest::Approx(1));

  Point3D ret(
    R_a*p.x+R_b*p.y+R_c*p.z,
    R_d*p.x+R_e*p.y+R_f*p.z,
    R_g*p.x+R_h*p.y+R_i*p.z
  );

  CHECK(ret.mag2()==doctest::Approx(1));

  return ret;
}

//Rotate the crap out of things
TEST_CASE("rotate"){
  for(int t=0;t<10;t++){
    std::vector<Point3D> pts;
    for(int i=0;i<100;i++)
      pts.push_back(Point2D(180.0-360.0*(rand()/(double)RAND_MAX),90.0-180.0*(rand()/(double)RAND_MAX)).toRadians().toXYZ(1));

    //Ensure everything is normalized before we continue
    for(const auto &p: pts)
      CHECK(p.mag2()==doctest::Approx(1));

    std::vector<double> angles;
    for(unsigned int i=0;i<pts.size();i++)
    for(unsigned int j=i+1;j<pts.size();j++)
      angles.emplace_back(pts[i].dot(pts[j]));

    for(int rn=0;rn<100;rn++){
      Point2D oldv(180.0-360.0*(rand()/(double)RAND_MAX),90.0-180.0*(rand()/(double)RAND_MAX));
      Point2D newv(180.0-360.0*(rand()/(double)RAND_MAX),90.0-180.0*(rand()/(double)RAND_MAX));
      Rotator r(oldv.toRadians().toXYZ(1),newv.toRadians().toXYZ(1));
      for(auto &p: pts)
        p = r(p);
      int count = 0;
      for(unsigned int i=0;i<pts.size();i++)
      for(unsigned int j=i+1;j<pts.size();j++)
        CHECK(angles[count++]==doctest::Approx(pts[i].dot(pts[j])));
    }
  }
}



TEST_CASE("Point2D: Constructor"){
  Point2D p;
}

Point2D::Point2D(double x0, double y0) {
  x = x0;
  y = y0;
}



Point2D& Point2D::toRadians() {
  assert(-180<=x && x<=180);
  assert(-90<=y && y<=90);

  x *= DEG_TO_RAD;
  y *= DEG_TO_RAD;

  assert(-M_PI<=x && x<=M_PI);
  assert(-M_PI/2<=y && y<=M_PI/2);
  return *this;
}

TEST_CASE("Point2D: Conversion to Radians"){
  Point2D p(-93, 45);
  p.toRadians();
  CHECK(p.x==doctest::Approx(-93*DEG_TO_RAD));
  CHECK(p.y==doctest::Approx(45*DEG_TO_RAD));
}



Point2D& Point2D::toDegrees() {
  assert(-M_PI<=x && x<=M_PI);
  assert(-M_PI/2<=y && y<=M_PI/2);

  x *= RAD_TO_DEG;
  y *= RAD_TO_DEG;

  assert(-180<=x && x<=180);
  assert(-90<=y && y<=90);
  return *this;
}

TEST_CASE("Point2D: Conversion to Degrees"){
  Point2D p(-93*DEG_TO_RAD, 45*DEG_TO_RAD);
  p.toDegrees();
  CHECK(p.x==doctest::Approx(-93));
  CHECK(p.y==doctest::Approx(45));
}



Point2D& Point2D::rotateTheta(const double rtheta) {
  x += M_PI+rtheta;        //Move [0,360] and add theta
  x = std::fmod(x,2*M_PI); //Map back to [0,360]
  if(x<0)
    x += 2*M_PI;
  x -= M_PI;               //Move back to [-180,180] system
  return *this;
}



bool Point2D::operator==(const Point2D &other) const {
  return x==other.x && y==other.y;
}



TEST_CASE("rotateTheta"){
  SUBCASE("Test1"){
    auto p = Point2D(-93,45).toRadians();
    p.rotateTheta(93*DEG_TO_RAD);
    CHECK(p.x==doctest::Approx(0));
  }
  SUBCASE("Test2"){
    auto p = Point2D(0,90).toRadians();
    p.rotateTheta(-67*DEG_TO_RAD);
    CHECK(p.x==doctest::Approx(-67*DEG_TO_RAD));
  }
  SUBCASE("Test3"){
    auto p = Point2D(144,90).toRadians();
    p.rotateTheta(72*DEG_TO_RAD);
    CHECK(p.x==doctest::Approx(-144*DEG_TO_RAD));
  }
}



Point3D Point2D::toXYZ(const double radius) const {
  assert(-M_PI<=x && x<=M_PI);
  assert(-M_PI/2<=y && y<=M_PI/2);

  const auto ys = M_PI-(y+M_PI/2); //Shift from [-90,90] -> [0,180]. 0 is North in the transformed coordinates.

  auto ret = Point3D(
    radius * std::cos(x) * std::sin(ys),
    radius * std::sin(x) * std::sin(ys),
    radius * std::cos(ys)
  );

  assert(ret.mag2());

  return ret;
}

TEST_CASE("Point 2D: toXYZ"){
  SUBCASE("North Pole"){
    const auto p3 = Point2D(0, 90).toRadians().toXYZ(1);
    CHECK(p3.x==doctest::Approx(0));
    CHECK(p3.y==doctest::Approx(0));
    CHECK(p3.z==doctest::Approx(1));
  }

  SUBCASE("South Pole"){
    const auto p3 = Point2D(0, -90).toRadians().toXYZ(1);
    CHECK(p3.x==doctest::Approx(0));
    CHECK(p3.y==doctest::Approx(0));
    CHECK(p3.z==doctest::Approx(-1));
  }
}



Point3D::Point3D(double x0, double y0, double z0) {
  x = x0;
  y = y0;
  z = z0;
}

Point2D Point3D::toLatLon() const {
  const double radius = mag();

  if(!(std::abs(radius-1)<1e-6 || std::abs(radius-6371)<1e-6))
    throw std::runtime_error("Radius invalid = " + std::to_string(radius));

  double x2 = std::atan2(y,x);
  double y2 = std::acos(z/radius);

  y2  = M_PI/2-y2; //From 0-up [0,180] to 90 up [-90,90]

  return Point2D(x2,y2);
}

TEST_CASE("Point 3D: toLatLon"){
  SUBCASE("North Pole"){
    const auto p2 = Point3D(0,0,1).toLatLon();
    //p2.x is not uniquely defined in the context of this test, so leave it out
    CHECK(p2.y==doctest::Approx(90*DEG_TO_RAD));
  }

  SUBCASE("North Pole Lat Lon"){
    const auto p2 = Point2D(0,90).toRadians().toXYZ(1).toLatLon();
    //p2.x is not uniquely defined in the context of this test, so leave it out
    CHECK(p2.y==doctest::Approx(90*DEG_TO_RAD));
  }

  SUBCASE("South Pole"){
    const auto p2 = Point3D(0, 0, -1).toLatLon();
    //p2.x is not uniquely defined in the context of this test, so leave it out
    CHECK(p2.y==doctest::Approx(-90*DEG_TO_RAD));
  }

  SUBCASE("South Pole Lat Lon"){
    const auto p2 = Point2D(0,-90).toRadians().toXYZ(1).toLatLon();
    //p2.x is not uniquely defined in the context of this test, so leave it out
    CHECK(p2.y==doctest::Approx(-90*DEG_TO_RAD));
  }

  SUBCASE("East Pole"){
    const auto p2 = Point2D(0, 0).toRadians().toXYZ(1).toLatLon();
    CHECK(p2.x==doctest::Approx(0));
    CHECK(p2.y==doctest::Approx(0));
  }

  SUBCASE("East Pole with 6371 radius"){
    const auto p2 = Point2D(0, 0).toRadians().toXYZ(6371).toLatLon();
    CHECK(p2.x==doctest::Approx(0));
    CHECK(p2.y==doctest::Approx(0));
  }

  SUBCASE("Minneapolis LatLon->XYZ->LatLon"){
    const auto p3 = Point2D(-93, 45).toRadians().toXYZ(1).toLatLon().toDegrees();
    CHECK(p3.x==doctest::Approx(-93).epsilon(0.005));
    CHECK(p3.y==doctest::Approx(45).epsilon(0.005));
  }

  SUBCASE("Check cube"){
    double lons[4] = {-135,-45,45,135};
    double lats[2] = {45,-45};
    for(int y=0;y<2;y++)
    for(int x=0;x<4;x++){
      auto p2 = Point2D(lons[x],lats[y]).toRadians().toXYZ(1).toLatLon().toDegrees();
      CHECK(lons[x]==doctest::Approx(p2.x));
      CHECK(lats[y]==doctest::Approx(p2.y));
    }
  }

  SUBCASE("Large numbers of random points with R=1"){
    for(int i=0;i<1000;i++){
      const double lon = 180.0-360.0*(rand()/(double)RAND_MAX);
      const double lat = 90.0-180.0*(rand()/(double)RAND_MAX);
      assert(-180<=lon && lon<=180);
      assert(-90<=lat && lat<=90);
      auto p2 = Point2D(lon,lat).toRadians().toXYZ(1).toLatLon().toDegrees();
      CHECK(lon==doctest::Approx(p2.x));
      CHECK(lat==doctest::Approx(p2.y));
    }
  }

  SUBCASE("Large numbers of random points with R=6371"){
    for(int i=0;i<1000;i++){
      const double lon = 180.0-360.0*(rand()/(double)RAND_MAX);
      const double lat = 90.0-180.0*(rand()/(double)RAND_MAX);
      assert(-180<=lon && lon<=180);
      assert(-90<=lat && lat<=90);
      const auto p2 = Point2D(lon,lat).toRadians().toXYZ(6371).toLatLon().toDegrees();
      CHECK(lon==doctest::Approx(p2.x));
      CHECK(lat==doctest::Approx(p2.y));
    }
  }
}



double Point3D::dot(const Point3D &b) const {
  return x*b.x+y*b.y+z*b.z;
}

Point3D Point3D::cross(const Point3D &b) const {
  Point3D temp;
  temp.x = -y*b.z + z*b.y;
  temp.y =  x*b.z - z*b.x;
  temp.z = -x*b.y + y*b.x;
  return temp;
}

Point3D& Point3D::unitify() {
  const double len = mag();
  if(len==0)
    return *this;
  x /= len;
  y /= len;
  z /= len;
  return *this;
}

double Point3D::mag2() const {
  return x*x+y*y+z*z;
}

double Point3D::mag() const {
  return std::sqrt(x*x+y*y+z*z);
}



TEST_CASE("3D Math"){
  CHECK(Point3D(1,0,0).dot(Point3D(0,1,0))==0);
  CHECK(Point3D(1,0,0).dot(Point3D(1,0,0))==1);
}
