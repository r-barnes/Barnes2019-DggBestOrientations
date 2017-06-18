#include "POI.hpp"
#include "doctest.h"

Orientation::Orientation(const Point2D &pole0, const double theta0){
  pole  = pole0;
  theta = theta0;
}

OrientationWithStats::OrientationWithStats(const Orientation &o) : Orientation(o) {}

TEST_CASE("Orientation Constructors"){
  Orientation o(Point2D(-93,45),13.7);
  Orientation ob(o);
  OrientationWithStats ows(o);
  CHECK(o.pole.x==ob.pole.x);
  CHECK(o.theta==ob.theta);
  CHECK(o.pole.x==ows.pole.x);
  CHECK(o.theta==ows.theta);
}