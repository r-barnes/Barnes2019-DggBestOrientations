#include "Point.hpp"
#include <cmath>
#include <cassert>
#include <cstdlib>
#include "doctest.h"

const double DEG_TO_RAD = M_PI/180.0;
const double RAD_TO_DEG = 180.0/M_PI;

Rotator::Rotator(const Point3D &oldv, const Point3D &newv){
  //Ensure that incoming vectors are normalized
  CHECK(oldv.x*oldv.x+oldv.y*oldv.y+oldv.z*oldv.z==doctest::Approx(1));
  CHECK(newv.x*newv.x+newv.y*newv.y+newv.z*newv.z==doctest::Approx(1));

  //Icosahedron's North pole
  const double a1 = oldv.x;
  const double a2 = oldv.y;
  const double a3 = oldv.z;

  //Vector defining the new pole
  const double b1 = newv.x;
  const double b2 = newv.y;
  const double b3 = newv.z;

  //B x A: This vector is perpendicular to both the North pole and the new pole
  double r1 = -a2*b3 + a3*b2;
  double r2 =  a1*b3 - a3*b1;
  double r3 = -a1*b2 + a2*b1;
  mr        = std::sqrt(r1*r1+r2*r2+r3*r3); //mr used in operator()

  //If the vector cannot be normalized it means the North Pole was coincident
  //with the new pole
  if(mr==0)
    return;

  //Normalize the vector
  r1 /= mr;
  r2 /= mr;
  r3 /= mr;

  //Ensure that things were actually normalized
  CHECK(r1*r1+r2*r2+r3*r3==doctest::Approx(1));

  const double c = a1*b1 + a2*b2 + a3*b3;  //cos theta = B dot A
  const double s = std::sqrt(1-c*c);            //= sin theta
  const double t = 1-c;

  CHECK(s*s+c*c==doctest::Approx(1));

  //Rotation Matrix from
  //Glassner, A.S. (Ed.), 1993. Graphics Gems I, 1st ed. p. 466
  R_a = c + r1*r1*t;
  R_b = r1*r2*t + r3*s;
  R_c = r1*r3*t - r2*s;
  R_d = r1*r2*t - r3*s;
  R_e = c + r2*r2*t;
  R_f = r1*s + r2*r3*t;
  R_g = r1*r3*t + r2*s;
  R_h = -r1*s + r2*r3*t;
  R_i = c + r3*r3*t;
}

Point3D Rotator::operator()(const Point3D &p) const {
  if(mr==0)
    return p;

  CHECK(std::sqrt(p.x*p.x+p.y*p.y+p.z*p.z)==doctest::Approx(1));

  Point3D ret(
    R_a*p.x+R_b*p.y+R_c*p.z,
    R_d*p.x+R_e*p.y+R_f*p.z,
    R_g*p.x+R_h*p.y+R_i*p.z
  );

  CHECK(std::sqrt(ret.x*ret.x+ret.y*ret.y+ret.z*ret.z)==doctest::Approx(1));

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
      CHECK(std::sqrt(p.x*p.x+p.y*p.y+p.z*p.z)==doctest::Approx(1));

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
  x += M_PI+rtheta;                                  //Move [0,360] and add theta
  x  = std::fmod(2*M_PI+std::fmod(x,2*M_PI),2*M_PI); //Map back to [0,360]
  x -= M_PI;                                         //Move back to [-180,180] system
  return *this;
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

  assert(ret.x*ret.x+ret.y*ret.y+ret.z*ret.z);

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
  const double radius = std::sqrt(x*x + y*y + z*z);
  
  assert(std::abs(radius-1)<1e-6 || std::abs(radius-6371)<1e-6);

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

TEST_CASE("Point3D::dot"){
  CHECK(Point3D(1,0,0).dot(Point3D(0,1,0))==0);
  CHECK(Point3D(1,0,0).dot(Point3D(1,0,0))==1);
}
