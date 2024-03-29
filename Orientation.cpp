#include "Orientation.hpp"
#include "doctest.h"

#include <cmath>
#include <iostream>
#include <utility>

const double DEG_TO_RAD = M_PI/180.0;

Orientation::Orientation(const Point2D &pole0, const double theta0){
  pole  = pole0;
  theta = theta0;
}

OrientationWithStats::OrientationWithStats(const Orientation &o) : Orientation(o) {}

OrientationWithStats::OrientationWithStats(const Point2D &pole0, const double theta0) : Orientation(pole0, theta0) {}

TEST_CASE("Orientation Constructors"){
  Orientation o(Point2D(-93,45),13.7);
  Orientation ob(o);
  OrientationWithStats ows(o);
  CHECK(o.pole.x==ob.pole.x);
  CHECK(o.theta==ob.theta);
  CHECK(o.pole.x==ows.pole.x);
  CHECK(o.theta==ows.theta);
}



OrientationGenerator::OrientationGenerator(
  const double point_spacingkm,
  const double south_of0,
  const double north_of0
) : south_of(south_of0*M_PI/180), north_of(north_of0*M_PI/180) {
  //Number of points to sample
  N = (long)(8*M_PI*Rearth*Rearth/std::sqrt(3)/point_spacingkm/point_spacingkm);

  // #pragma omp critical
  // {
  //   std::cerr << "Initializing orientation generator"      <<std::endl;
  //   std::cerr << "\tpoint_spacingkm = " << point_spacingkm <<std::endl;
  //   std::cerr << "\tradial_limit    = " << radial_limit    <<std::endl;
  //   std::cerr << "\tN               = " << N               <<std::endl;
  //   std::cerr << "\tNmax            = " << Nmax            <<std::endl;
  // }
}



///Returns the number of orientations to be generated
long OrientationGenerator::size() const {
  return N;
}

//Generate a set of orientations. The zeroth pole is near the South Pole with
//orientations spiraling outwards from there until the `radial_limit` is
//reached: the maximum number of radians outward from the South Pole. Poles are
//spaced with approximately `point_spacingkm` kilometres between themselves. The
//polyhedron will also be rotated from `theta_min` to `theta_max` with steps of
//size `theta_step`.
std::pair<bool, Point2D> OrientationGenerator::operator()(long i) const {
  Point2D pole (
    M_PI*(3.0-std::sqrt(5.0))*i,
    std::acos(1-(2.0*i+1.0)/N)
  );
  pole.x = std::fmod(pole.x,2*M_PI)-M_PI;
  pole.y = pole.y-M_PI/2;
  pole.y = -pole.y; //Orientate so North Pole is up

  const bool valid = north_of <= pole.y && pole.y <= south_of;

  return std::make_pair(valid, pole);
}
