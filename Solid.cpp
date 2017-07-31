#include "Solid.hpp"
#include "GeoStuff.hpp"
#include <iomanip>
#include <iostream>
#include <limits>
#include <cassert>
#include "doctest.h"

const double DEG_TO_RAD = M_PI/180.0;
const double RAD_TO_DEG = 180.0/M_PI;

const double phi = (1+std::sqrt(5))/2;



SolidXY& SolidXY::rotate(const Point2D &p, double rtheta){
  rotateTheta(rtheta);
  *this = toXYZ(1).rotateTo(p.toXYZ(1)).toLatLon();
  return *this;
}



SolidXY& SolidXY::rotate(double rlat, double rlon, double rtheta){
  Point2D temp(rlon,rlat);
  return rotate(temp,rtheta);
}



SolidXY& SolidXY::rotateTheta(const double rtheta){
  if(rtheta==0)
    return *this;
  
  assert(std::abs(v[0].y-M_PI/2)<1e-6); //Can only rotate when icosahedron is North-South aligned

  //0 and 1 are the poles, which do not rotate
  for(unsigned int i=2;i<v.size();i++)
    v[i].rotateTheta(rtheta);

  return *this;
}



SolidXY& SolidXY::toMercator(){
  for(auto &p: v)
    p = WGS84toEPSG3857(p);
  return *this;
}



SolidXY& SolidXY::toRadians(){
  for(auto &p: v)
    p.toRadians();
  return *this;
}



//Return neighbours in a N1a,N1b,N2a,N2b,... format
std::vector<int> SolidXY::neighbors() const {
  std::vector<int> ret;
  double dist = std::numeric_limits<double>::infinity();
  //Find nearest vertex that isn't itself
  for(unsigned int i=1;i<v.size();i++)
    dist = std::min(dist,GeoDistanceHaversine(v[0],v[i]));
  //Increase by 5% to account for floating point errors
  dist *= 1.3;
  for(unsigned int i=0;  i<v.size();i++)
  for(unsigned int j=i+1;j<v.size();j++)
    if(GeoDistanceHaversine(v[i],v[j])<dist){ //SolidXY is about as close as any close pole
      ret.emplace_back(i);
      ret.emplace_back(j);
    }
  return ret;
}



double SolidXY::neighborDistance() const {
  const auto n = neighbors();
  return GeoDistanceHaversine(v[n[0]], v[n[1]]);
}



SolidXYZ SolidXY::toXYZ(const double radius) const {
  SolidXYZ temp;
  temp.v.resize(v.size());
  for(unsigned int i=0;i<v.size();i++)
    temp.v[i] = v[i].toXYZ(radius);
  return temp;
}



SolidXY SolidXYZ::toLatLon() const {
  SolidXY temp;
  temp.v.resize(v.size());
  for(unsigned int i=0;i<v.size();i++)
    temp.v[i] = v[i].toLatLon();
  return temp;
}



//https://math.stackexchange.com/a/476311/14493
SolidXYZ& SolidXYZ::rotateTo(const Point3D &o){
  Rotator r(v[0],o);

  rotate(r);

  return *this;
}



SolidXYZ& SolidXYZ::normalize() {
  const double dist_to_origin = std::sqrt(v[0].x*v[0].x + v[0].y*v[0].y + v[0].z*v[0].z);
  for(auto &p: v){
    p.x /= dist_to_origin;
    p.y /= dist_to_origin;
    p.z /= dist_to_origin;
  }
  return *this;
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



SolidXYZ& SolidXYZ::rotate(const Rotator &r) {
  //Rotate each pole
  for(auto &p: v)
    p = r(p);

  return *this;
}



SolidXY OrientationToRegularIcosahedron(const Orientation &o) {
  const double IEL = std::atan(0.5); //Icosahedron equatorial latitudes
  const double IES = 36*DEG_TO_RAD;    //Icosahedron equatorial spacing

  SolidXY sxy;
  sxy.v = {{
    { 0    ,  M_PI/2},
    { 0    , -M_PI/2},
    {-5*IES,     IEL},
    {-4*IES,    -IEL},
    {-3*IES,     IEL},
    {-2*IES,    -IEL},
    {-1*IES,     IEL},
    { 0*IES,    -IEL},
    { 1*IES,     IEL},
    { 2*IES,    -IEL},
    { 3*IES,     IEL},
    { 4*IES,    -IEL}
  }};

  sxy.rotate(o.pole,o.theta);
  return sxy;
}


SolidXY OrientationToRegularDodecahedron(const Orientation &o){
  static std::vector<Point2D> v;
  if(v.size()==0){
    const double phii = 1/phi;
    SolidXYZ sxyz;
    sxyz.v = {
      { 1, 1, 1},
      { 1, 1,-1},
      { 1,-1, 1},
      { 1,-1,-1},
      {-1, 1, 1},
      {-1, 1,-1},
      {-1,-1, 1},
      {-1,-1,-1},
      { 0, phi, phii},
      { 0, phi,-phii},
      { 0,-phi, phii},
      { 0,-phi,-phii},
      { phii,0, phi},
      { phii,0,-phi},
      {-phii,0, phi},
      {-phii,0,-phi},
      { phi, phii,0},
      { phi,-phii,0},
      {-phi, phii,0},
      {-phi,-phii,0}
    };
    sxyz.normalize();
    Rotator r(sxyz.v[0],Point3D(0,0,1));
    sxyz.rotate(r);
    auto sxy = sxyz.toLatLon();
    v = sxy.v;
  }

  SolidXY sxy;
  sxy.v = v;
  sxy.rotate(o.pole,o.theta);

  return sxy;
}

SolidXY OrientationToRegularTetrahedron(const Orientation &o){
  static std::vector<Point2D> v;
  if(v.size()==0){
    SolidXYZ sxyz;
    sxyz.v = {
      { 1, 0,-1/std::sqrt(2)},
      {-1, 0,-1/std::sqrt(2)},
      { 0, 1, 1/std::sqrt(2)},
      { 0,-1, 1/std::sqrt(2)}
    };
    sxyz.normalize();
    Rotator r(sxyz.v[0],Point3D(0,0,1));
    sxyz.rotate(r);
    auto sxy = sxyz.toLatLon();
    v        = sxy.v;
  }

  SolidXY sxy;
  sxy.v = v;
  sxy.rotate(o.pole,o.theta);

  return sxy;
}

SolidXY OrientationToRegularOctahedron(const Orientation &o){
  static std::vector<Point2D> v;
  if(v.size()==0){
    SolidXYZ sxyz;
    sxyz.v = {
      { 1,0,0},
      {-1,0,0},
      {0, 1,0},
      {0,-1,0},
      {0,0, 1},
      {0,0,-1}
    };
    sxyz.normalize();
    Rotator r(sxyz.v[0],Point3D(0,0,1));
    sxyz.rotate(r);
    auto sxy = sxyz.toLatLon();
    v        = sxy.v;
  }

  SolidXY sxy;
  sxy.v = v;
  sxy.rotate(o.pole,o.theta);
  
  return sxy;
}

SolidXY OrientationToCuboctahedron(const Orientation &o){
  static std::vector<Point2D> v;
  if(v.size()==0){
    SolidXYZ sxyz;
    sxyz.v = {
      { 1, 1, 0},
      { 1,-1, 0},
      {-1, 1, 0},
      {-1,-1, 0},
      { 1, 0, 1},
      { 1, 0,-1},
      {-1, 0, 1},
      {-1, 0,-1},
      { 0, 1, 1},
      { 0, 1,-1},
      { 0,-1, 1},
      { 0,-1,-1}
    };
    sxyz.normalize();
    Rotator r(sxyz.v[0],Point3D(0,0,1));
    sxyz.rotate(r);
    auto sxy = sxyz.toLatLon();
    v        = sxy.v;
  }

  SolidXY sxy;
  sxy.v = v;
  sxy.rotate(o.pole,o.theta);
  
  return sxy;
}

SolidXY OrientationFullerIcosahedron(){
  static std::vector<Point2D> v;
  if(v.size()==0){
    v = {{
      {  10.53620,  64.7     },
      {  -5.24539,   2.300882},
      {  58.15771,  10.447378},
      { 122.3    ,  39.1     },
      {-143.47849,  50.103201},
      { -67.13233,  23.717925},
      { -57.7    , -39.1     },
      {  36.5215 , -50.1032  },
      { 112.86767, -23.717925},
      { 174.7546 ,  -2.3009  },
      {-121.84229, -10.447345},
      {-169.4638 , -64.7     }
    }};
    SolidXY sxy;
    sxy.v = v;
    sxy.toRadians();
    v = sxy.v;
  }

  SolidXY sxy;
  sxy.v = v;

  return sxy;
}

SolidXY OrientationToPoint(const Orientation &o){
  SolidXY sxy;
  sxy.v = {{o.pole}};
  return sxy;
}






TEST_CASE("Shape metrics"){
  const Orientation o(Point2D(0,90).toRadians(),0);

  auto icosahedron         = OrientationToRegularIcosahedron(o);
  auto regulardodecahedron = OrientationToRegularDodecahedron(o);
  auto regulartetrahedron  = OrientationToRegularTetrahedron(o);
  auto regularoctahedron   = OrientationToRegularOctahedron(o);
  auto cuboctahedron       = OrientationToCuboctahedron(o);
  auto point               = OrientationToPoint(o);

  SUBCASE("Number of vertices"){
    SUBCASE("icosahedron")         {CHECK(icosahedron.v.size()        ==12 );}
    SUBCASE("regulardodecahedron") {CHECK(regulardodecahedron.v.size()==20 );}
    SUBCASE("regulartetrahedron")  {CHECK(regulartetrahedron.v.size() == 4 );}
    SUBCASE("regularoctahedron")   {CHECK(regularoctahedron.v.size()  == 6 );}
    SUBCASE("cuboctahedron")       {CHECK(cuboctahedron.v.size()      ==12 );}
    SUBCASE("point")               {CHECK(point.v.size()              == 1 );}
  } 

  SUBCASE("Equal Distances from Center"){
    const auto dist = [](const Point3D &p3d){return std::sqrt(p3d.x*p3d.x+p3d.y*p3d.y+p3d.z*p3d.z);};
    const auto all_dist_equal = [&](const SolidXYZ &sxyz){
      const auto dist_nominal = dist(sxyz.v[0]);
      for(const auto &p: sxyz.v)
        CHECK(dist(p)==doctest::Approx(dist_nominal));
    };
    SUBCASE("icosahedron")         {all_dist_equal(icosahedron.toXYZ(1));        }
    SUBCASE("regulardodecahedron") {all_dist_equal(regulardodecahedron.toXYZ(1));}
    SUBCASE("regulartetrahedron")  {all_dist_equal(regulartetrahedron.toXYZ(1)); }
    SUBCASE("regularoctahedron")   {all_dist_equal(regularoctahedron.toXYZ(1));  }
    SUBCASE("cuboctahedron")       {all_dist_equal(cuboctahedron.toXYZ(1));      }
    SUBCASE("point")               {all_dist_equal(point.toXYZ(1));              }
  }

  SUBCASE("Point at North Pole"){
    SUBCASE("icosahedron"){
      CHECK(icosahedron.v[0].x==doctest::Approx(0));
      CHECK(icosahedron.v[0].y==doctest::Approx(M_PI/2));
    }
    SUBCASE("regulardodecahedron"){
      CHECK(regulardodecahedron.v[0].x==doctest::Approx(0));
      CHECK(regulardodecahedron.v[0].y==doctest::Approx(M_PI/2));
    }
    SUBCASE("regulartetrahedron"){
      CHECK(regulartetrahedron.v[0].x==doctest::Approx(0));
      CHECK(regulartetrahedron.v[0].y==doctest::Approx(M_PI/2));
    }
    SUBCASE("regularoctahedron"){
      CHECK(regularoctahedron.v[0].x==doctest::Approx(0));
      CHECK(regularoctahedron.v[0].y==doctest::Approx(M_PI/2));
    }
    SUBCASE("cuboctahedron"){
      CHECK(cuboctahedron.v[0].x==doctest::Approx(0));
      CHECK(cuboctahedron.v[0].y==doctest::Approx(M_PI/2));
    }
    SUBCASE("point"){
      CHECK(point.v[0].x==doctest::Approx(0));
      CHECK(point.v[0].y==doctest::Approx(M_PI/2));
    }
  }

  SUBCASE("Equidistant"){
    const auto dist_checker = [](const SolidXY &sxy){
      const auto neighbors = sxy.neighbors();
      std::vector<int> ncount(sxy.v.size(),0);
      for(const auto n: neighbors)
        ncount[n]++;
      const int ncount_nominal = ncount[0];
      for(const auto nc: ncount)
        CHECK(nc==ncount_nominal);
    };
    SUBCASE("icosahedron")         {dist_checker(icosahedron);         }
    SUBCASE("regulardodecahedron") {dist_checker(regulardodecahedron); }
    SUBCASE("regulartetrahedron")  {dist_checker(regulartetrahedron);  }
    SUBCASE("regularoctahedron")   {dist_checker(regularoctahedron);   }
    SUBCASE("cuboctahedron")       {dist_checker(cuboctahedron);       }
    SUBCASE("point")               {dist_checker(point);               }
  }

  SUBCASE("Edge counts"){
    const auto count_edges = [](const SolidXY &sxy){
      return sxy.toXYZ(1).neighbors().size()/2; //Two neighbours equals one edge
    };

    SUBCASE("icosahedron")         {CHECK(count_edges(icosahedron)         == 30 );}
    SUBCASE("regulardodecahedron") {CHECK(count_edges(regulardodecahedron) == 30 );}
    SUBCASE("regulartetrahedron")  {CHECK(count_edges(regulartetrahedron)  ==  6 );}
    SUBCASE("regularoctahedron")   {CHECK(count_edges(regularoctahedron)   == 12 );}
    SUBCASE("cuboctahedron")       {CHECK(count_edges(cuboctahedron)       == 24 );}
    SUBCASE("point")               {CHECK(count_edges(point)               ==  0 );}
  }

  SUBCASE("Neighbour counts and distances"){
    const auto ncounter = [](const SolidXY &sxy){
      const auto n = sxy.toXYZ(1).neighbors();

      const auto ndist = sxy.neighborDistance();

      std::vector< std::vector<int> > ncheck(sxy.v.size());
      for(unsigned int ni=0;ni<n.size();ni+=2){
        ncheck.at(n.at(ni)).emplace_back(n.at(ni+1));
        ncheck.at(n.at(ni+1)).emplace_back(n.at(ni));
        CHECK(GeoDistanceHaversine(sxy.v.at(n.at(ni)),sxy.v.at(n.at(ni+1)))==doctest::Approx(ndist));
      }

      //Choose one vertex randomly. All should match it.
      const auto n_to_match = ncheck.front().size();

      //Each vertex should have 5 neighbours
      for(const auto &ni: ncheck)
        CHECK(ni.size()==n_to_match);

      return n_to_match;
    };

    SUBCASE("icosahedron")         {CHECK(ncounter(icosahedron)         == 5 );}
    SUBCASE("regulardodecahedron") {CHECK(ncounter(regulardodecahedron) == 3 );}
    SUBCASE("regulartetrahedron")  {CHECK(ncounter(regulartetrahedron)  == 3 );}
    SUBCASE("regularoctahedron")   {CHECK(ncounter(regularoctahedron)   == 4 );}
    SUBCASE("cuboctahedron")       {CHECK(ncounter(cuboctahedron)       == 4 );}
  }
}

TEST_CASE("rotate"){
  const Orientation o(Point2D(0,90).toRadians(),0);
  const auto a = OrientationToRegularIcosahedron(o);
  const auto p = OrientationToRegularIcosahedron(o).rotate(Point2D(0,90).toRadians(),0);
  for(unsigned int i=0;i<a.v.size();i++){
    CHECK(a.v[i].x==doctest::Approx(p.v[i].x));
    CHECK(a.v[i].y==doctest::Approx(p.v[i].y));
  }
}


TEST_CASE("rotateTheta"){
  const Orientation o(Point2D(0,90).toRadians(),0);
  SolidXY s = OrientationToRegularIcosahedron(o);
  CHECK(s.v.at(0).x==doctest::Approx(0));
  CHECK(s.v.at(7).x==doctest::Approx(0));
  CHECK(s.v.at(11).x==doctest::Approx(4*36*DEG_TO_RAD));
  s.rotateTheta(72*DEG_TO_RAD);
  CHECK(s.v.at(0).x==doctest::Approx(0));
  CHECK(s.v.at(7).x==doctest::Approx(72*DEG_TO_RAD));
  CHECK(s.v.at(11).x==doctest::Approx(-144*DEG_TO_RAD));
}


TEST_CASE("toMercator"){
  const Orientation o(Point2D(0,90).toRadians(),0);
  SolidXY ico2d = OrientationToRegularIcosahedron(o);
  auto p     = WGS84toEPSG3857(ico2d.v[3]);
  ico2d.toMercator();
  CHECK(p.x==ico2d.v[3].x);
  CHECK(p.y==ico2d.v[3].y);
}

TEST_CASE("toRadians"){
  const Orientation o(Point2D(0,90).toRadians(),0);
  SolidXY ico2d = OrientationToRegularIcosahedron(o);
  auto p     = ico2d.v[3];
  ico2d.toRadians();
  p.toRadians();
  CHECK(p.x==ico2d.v[3].x);
  CHECK(p.y==ico2d.v[3].y);
}

TEST_CASE("toXYZ"){
  const Orientation o(Point2D(0,90).toRadians(),0);
  SolidXY a = OrientationToRegularIcosahedron(o);
  SolidXY b = OrientationToRegularIcosahedron(o).toXYZ(1).toLatLon();
  for(unsigned int i=0;i<a.v.size();i++){
    CHECK(a.v[i].x==doctest::Approx(b.v[i].x));
    CHECK(a.v[i].y==doctest::Approx(b.v[i].y));
  }
}

TEST_CASE("SolidXYZ::toLatLon"){
  const Orientation o(Point2D(0,90).toRadians(),0);
  SolidXY a = OrientationToRegularIcosahedron(o);
  SolidXY b = OrientationToRegularIcosahedron(o).toXYZ(1).toLatLon();
  for(unsigned int i=0;i<a.v.size();i++){
    CHECK(a.v[i].x==doctest::Approx(b.v[i].x));
    CHECK(a.v[i].y==doctest::Approx(b.v[i].y));
  }
}

TEST_CASE("rotateTo"){
  SUBCASE("theta rotation"){
    const Orientation o(Point2D(0,90).toRadians(),0);
    SolidXY a = OrientationToRegularIcosahedron(o);
    SolidXY r = OrientationToRegularIcosahedron(o).rotate(90*DEG_TO_RAD,0,45*DEG_TO_RAD); //Rotate forward
    r.rotate(90*DEG_TO_RAD,0,-45*DEG_TO_RAD);                                      //Rotate back
    for(unsigned int i=0;i<a.v.size();i++){
      CHECK(a.v[i].x==doctest::Approx(r.v[i].x));
      CHECK(a.v[i].y==doctest::Approx(r.v[i].y));
    }
  }
}
