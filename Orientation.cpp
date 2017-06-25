#include "Orientation.hpp"
#include "doctest.h"
#include <iostream>
#include <cmath>

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
  const double radial_limit
){
  //Number of points to sample
  N = (long)(8*M_PI*Rearth*Rearth/std::sqrt(3)/point_spacingkm/point_spacingkm);

  //This might induce minor sample loss near the South Pole due to numeric
  //issues, but that turns out not to be an issue since we generate solids from
  //the North Pole and only explore orientations down to the equator (below
  //which symmetry guarantees that we've already explored what we need to)
  Nmax = (N*(1-std::cos(radial_limit))-1)/2;

  #pragma omp critical
  {
    std::cerr << "Initializing orientation generator"      <<std::endl;
    std::cerr << "\tpoint_spacingkm = " << point_spacingkm <<std::endl;
    std::cerr << "\tradial_limit    = " << radial_limit    <<std::endl;
    std::cerr << "\tN               = " << N               <<std::endl;
    std::cerr << "\tNmax            = " << Nmax            <<std::endl;
  }
}

long OrientationGenerator::size() const {
  return Nmax;
}

//Generate a set of orientations. The zeroth pole is near the South Pole with
//orientations spiraling outwards from there until the `radial_limit` is
//reached: the maximum number of radians outward from the South Pole. Poles are
//spaced with approximately `point_spacingkm` kilometres between themselves. The
//polyhedron will also be rotated from `theta_min` to `theta_max` with steps of
//size `theta_step`.
Point2D OrientationGenerator::operator()(long i) const {
  Point2D pole (
    M_PI*(3.0-std::sqrt(5.0))*i,
    std::acos(1-(2.0*i+1.0)/N)
  );
  pole.x = std::fmod(pole.x,2*M_PI)-M_PI;
  pole.y = pole.y-M_PI/2;
  pole.y = -pole.y; //Orientate so North Pole is up
  return pole;
}
