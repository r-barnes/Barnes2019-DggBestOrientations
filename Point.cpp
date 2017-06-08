#include "Point.hpp"
#include <cmath>
#include <cassert>
#include "doctest.h"

const double DEG_TO_RAD = M_PI/180.0;
const double RAD_TO_DEG = 180.0/M_PI;

Point2D::Point2D(){}

TEST_CASE("Point2D: Constructor"){
  Point2D p;
}

Point2D::Point2D(double x, double y) {
  this->x = x;
  this->y = y;
}



Point2D& Point2D::toRadians() {
  x *= DEG_TO_RAD;
  y *= DEG_TO_RAD;
  return *this;
}

TEST_CASE("Point2D: Conversion to Radians"){
  Point2D p(-93, 45);
  p.toRadians();
  CHECK(p.x==doctest::Approx(-93*DEG_TO_RAD));
  CHECK(p.y==doctest::Approx(45*DEG_TO_RAD));
}



Point2D& Point2D::toDegrees() {
  x *= RAD_TO_DEG;
  y *= RAD_TO_DEG;
  return *this;
}

TEST_CASE("Point2D: Conversion to Degrees"){
  Point2D p(-93*DEG_TO_RAD, 45*DEG_TO_RAD);
  p.toDegrees();
  CHECK(p.x==doctest::Approx(-93));
  CHECK(p.y==doctest::Approx(45));
}






Point3D Point2D::toXYZ(const double radius) const {
  auto ret = Point3D(
    radius * std::cos(x+M_PI) * std::sin(M_PI-(y+M_PI/2)),
    radius * std::sin(x+M_PI) * std::sin(M_PI-(y+M_PI/2)),
    radius * std::cos(M_PI-(y+M_PI/2))
  );
  assert(ret.x*ret.x+ret.y*ret.y+ret.z*ret.z);
  return ret;
}

TEST_CASE("Point 2D: toXYZ"){
  SUBCASE("North Pole"){
    const auto p3 = Point2D(0, 90*DEG_TO_RAD).toXYZ(1);
    CHECK(p3.x==doctest::Approx(0));
    CHECK(p3.y==doctest::Approx(0));
    CHECK(p3.z==doctest::Approx(1));
  }

  SUBCASE("South Pole"){
    const auto p3 = Point2D(0, -90*DEG_TO_RAD).toXYZ(1);
    CHECK(p3.x==doctest::Approx(0));
    CHECK(p3.y==doctest::Approx(0));
    CHECK(p3.z==doctest::Approx(-1));
  }

  SUBCASE("East Pole"){
    const auto p3 = Point2D(0, 0*DEG_TO_RAD).toXYZ(1);
    CHECK(p3.x==doctest::Approx(-1));
    CHECK(p3.y==doctest::Approx(0));
    CHECK(p3.z==doctest::Approx(0));
  }

  SUBCASE("East Pole with 6371 radius"){
    const auto p3 = Point2D(0, 0*DEG_TO_RAD).toXYZ(6371);
    CHECK(p3.x==doctest::Approx(-6371));
    CHECK(p3.y==doctest::Approx(0));
    CHECK(p3.z==doctest::Approx(0));
  }

  SUBCASE("East Pole with 6371 radius"){
    const auto p3 = Point2D(-93*DEG_TO_RAD, 45*DEG_TO_RAD).toXYZ(6371);
    CHECK(p3.x==doctest::Approx(235.7722950021));
    CHECK(p3.y==doctest::Approx(4498.8033881144));
    CHECK(p3.z==doctest::Approx(4504.9773029395));
  }
}



Point3D::Point3D(double x, double y, double z) {
  this->x = x;
  this->y = y;
  this->z = z;
}

Point2D Point3D::toLatLon() const {
  const double radius = x*x + y*y + z*z; //TODO: Square root?
  return Point2D(
    std::atan2(y,x)-M_PI,
    std::asin(z/radius)
  );
}

TEST_CASE("Point 3D: toLatLon"){
  SUBCASE("North Pole"){
    const auto p2 = Point2D(0, 90*DEG_TO_RAD).toXYZ(1).toLatLon();
    CHECK(p2.x==doctest::Approx(0));
    CHECK(p2.y==doctest::Approx(90*DEG_TO_RAD));
  }

  SUBCASE("South Pole"){
    const auto p2 = Point2D(0, -90*DEG_TO_RAD).toXYZ(1).toLatLon();
    CHECK(p2.x==doctest::Approx(0));
    CHECK(p2.y==doctest::Approx(-90*DEG_TO_RAD));
  }

  SUBCASE("East Pole"){
    const auto p2 = Point2D(0, 0*DEG_TO_RAD).toXYZ(1).toLatLon();
    CHECK(p2.x==doctest::Approx(0));
    CHECK(p2.y==doctest::Approx(0));
  }

  SUBCASE("East Pole with 6371 radius"){
    const auto p2 = Point2D(0, 0*DEG_TO_RAD).toXYZ(6371).toLatLon();
    CHECK(p2.x==doctest::Approx(0));
    CHECK(p2.y==doctest::Approx(0));
  }

  SUBCASE("East Pole with 6371 radius"){
    const auto p3 = Point2D(-93*DEG_TO_RAD, 45*DEG_TO_RAD).toXYZ(1).toLatLon();
    CHECK(p3.x==doctest::Approx(-93*DEG_TO_RAD).epsilon(0.005));
    CHECK(p3.y==doctest::Approx(45*DEG_TO_RAD).epsilon(0.005));
  }
}



