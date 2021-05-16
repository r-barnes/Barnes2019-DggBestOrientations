#ifndef DOCTEST_CONFIG_DISABLE
  #define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#endif

#include "CLI11.hpp"
#include "GeoStuff.hpp"
#include "Orientation.hpp"
#include "OrientationIndex.hpp"
#include "Point.hpp"
#include "PointCloud.hpp"
#include "Polygon.hpp"
#include "Progress.hpp"
#include "random.hpp"
#include "Solid.hpp"
#include "SpIndex.hpp"
#include "Timer.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "doctest.h"

const double Rearth = 6371; //km

const double DEG_TO_RAD = M_PI/180.0;
const double RAD_TO_DEG = 180.0/M_PI;

//Polyhedron and Projection configuration globals
SolidifyingFunc solidifier_func;
TransLLto3D_t TransLLto3D;
Trans3DtoLL_t Trans3DtoLL;
GeoDistance_t GeoDistance;
std::string chosen_polyhedron;
std::string chosen_projection;
GreatCircleGeneratorType great_circle_generator;

class PointWithStats {
 public:
  double      dist;
  Point2D     pt;
  std::string label;
  PointWithStats() = default;
  PointWithStats(double dist0, const Point2D &pt0, const std::string &label0){
    dist  = dist0;
    pt    = pt0;
    label = label0;
  }
};

void SetupForProjection(const std::string &projection_name){
  chosen_projection = projection_name;
  if(projection_name=="spherical"){
    TransLLto3D            = WGS84toSphericalCartesian;
    Trans3DtoLL            = SphericalCartesiantoWGS84;
    GeoDistance            = GeoDistanceSphere;
    great_circle_generator = GreatCircleGeneratorType::SIMPLE_SPHERICAL;
  } else if(projection_name=="ellipsoidal"){
    TransLLto3D            = WGS84toEllipsoidCartesian;
    Trans3DtoLL            = EllipsoidCartesiantoWGS84;
    GeoDistance            = GeoDistanceEllipsoid;
    great_circle_generator = GreatCircleGeneratorType::SIMPLE_SPHERICAL;
  } else if(projection_name=="haversine"){
    TransLLto3D            = WGS84toSphericalCartesian;
    Trans3DtoLL            = SphericalCartesiantoWGS84;
    GeoDistance            = GeoDistanceHaversine;
    great_circle_generator = GreatCircleGeneratorType::SIMPLE_SPHERICAL;
  } else {
    throw std::runtime_error("Unrecognized projection!");
  }
}



//Kind of an ugly function that sorts pts in a clockwise ordering based on a
//given center point
// void ClockwiseSorter(std::vector< std::pair<unsigned int, double> > &pts, const Point3D &center_xyz, const PointCloud &wgsp4c){
//   typedef std::pair< Point3D, std::pair<unsigned int, double> > p3did;

//   //Construct rotation matrix
//   const Rotator rotator(center_xyz, Point3D(0,0,1));

//   //Rotate all points so that they project nicely into the xy-plane
//   std::vector<p3did> pts3d;
//   for(const auto &pidx: pts)
//     pts3d.emplace_back(rotator(wgsp4c.pts.at(pidx.first)), pidx);

//   std::sort(pts3d.begin(), pts3d.end(), [&](const p3did &ai, const p3did &bi){
//     const auto &a = ai.first;
//     const auto &b = bi.first;

//     //Quick quadrant-based method
//     if (a.x - center_xyz.x >= 0 && b.x - center_xyz.x < 0)
//       return true;

//     if (a.x - center_xyz.x < 0 && b.x - center_xyz.x >= 0)
//       return false;

//     if (a.x - center_xyz.x == 0 && b.x - center_xyz.x == 0) {
//       if (a.y - center_xyz.y >= 0 || b.y - center_xyz.y >= 0)
//         return a.y > b.y;
//       return b.y > a.y;
//     }

//     //More expensive cross-product method
//     // compute the cross product of vectors (center_xyz -> a) x (center_xyz -> b)
//     int det = (a.x - center_xyz.x) * (b.y - center_xyz.y) - (b.x - center_xyz.x) * (a.y - center_xyz.y);
//     if (det < 0)
//       return true;
//     if (det > 0)
//       return false;

//     // points a and b are on the same line from the center_xyz
//     // check which point is closer to the center_xyz
//     return ai.second.second > bi.second.second;
//   });

//   pts.clear();
//   for(const auto &p: pts3d)
//     pts.push_back(p.second);
// }


PointCloud ReadPointCloudFromShapefile(std::string filename, std::string layer, const double interp_dist){
  std::cerr<<"Reading point cloud from shapefile..."<<std::endl;
  auto landmass_wgs84 = ReadShapefile(filename, layer);

  std::cerr<<"Converting WGS84 polygons to radians..."<<std::endl;
  for(auto &p: landmass_wgs84)
    p.toRadians();

  const auto poly_fill_point = [&](const Point2D &a, const Point2D &b){
    std::vector<Point2D> ret;
    // Do a quick flat-earth projection check to see if the edge is longer
    // than the allowed length. Since the distances are short this should
    // be safe.
    const auto quickdist = GeoDistanceFlatEarth(a,b);
    if(quickdist>interp_dist){
      const auto gcg = GreatArcFactory::make(great_circle_generator, a, b, interp_dist);
      for(unsigned int i=0;i<gcg->size();i++)
        ret.push_back((*gcg)(i));
    }
    ret.push_back(b);
    return ret;
  };

  const auto poly_filler = [&](const std::vector<Point2D> &poly){
    std::vector<Point2D> ret;
    for(unsigned int i=0;i<poly.size()-1;i++){
      auto pfp = poly_fill_point(poly[i],poly[i+1]);
      ret.insert(ret.end(),pfp.begin(),pfp.end());
    }
    auto pfp = poly_fill_point(poly.back(),poly.front());
    ret.insert(ret.end(),pfp.begin(),pfp.end());
    return ret;
  };

  std::cerr<<"Ensuring polygon edges are not too long..."<<std::endl;
  std::vector<Point2D> wgs84_ll_flat;
  #pragma omp declare reduction (merge : std::vector<Point2D> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
  #pragma omp parallel for default(none) schedule(static) shared(landmass_wgs84,poly_filler) reduction(merge: wgs84_ll_flat)
  for(unsigned int pi=0;pi<landmass_wgs84.size();pi++){
    const auto filled_poly = poly_filler(landmass_wgs84[pi].exterior);
    wgs84_ll_flat.insert(wgs84_ll_flat.end(),filled_poly.begin(),filled_poly.end());
  }

  landmass_wgs84.clear();
  landmass_wgs84.shrink_to_fit();

  std::cerr<<"Converting points to WGS84 Cartesian..."<<std::endl;
  std::vector<Point3D> wgs84_xyz;
  #pragma omp declare reduction (merge : std::vector<Point3D> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
  #pragma omp parallel for default(none) schedule(static) shared(TransLLto3D,wgs84_ll_flat) reduction(merge:wgs84_xyz)
  for(unsigned int i=0;i<wgs84_ll_flat.size();i++)
    wgs84_xyz.push_back(TransLLto3D(wgs84_ll_flat[i]));
    //wgs84_xyz.push_back(wgs84_ll_flat[i].toXYZ(1));

  wgs84_ll_flat.clear();
  wgs84_ll_flat.shrink_to_fit();

  wgs84_xyz.shrink_to_fit();

  std::cerr<<"Number of points = "<<wgs84_xyz.size()<<std::endl;

  PointCloud wgs84pc;
  std::cerr<<"Adding polygon points to PointCloud..."<<std::endl;
  wgs84pc.pts.swap(wgs84_xyz);

  return wgs84pc;
}



double DistanceToCoast(const Point2D &pt, const PointCloud &wgs84pc){
  const auto cp  = wgs84pc.queryPoint(TransLLto3D(pt)); //Closest point
  const auto llc = Trans3DtoLL(cp);
  return GeoDistance(llc,pt);
}



class HillClimber {
 public:
  static const int coords = 2;
  int steps               = 0;
  PointWithStats start;
  PointWithStats best;
 private:
  double wrapPolar(double x, const double rangemin, const double rangemax) const {
    if(x<rangemin)
      x = std::abs(x-rangemin)+rangemin;
    else if(x>rangemax)
      x = rangemax-std::abs(x-rangemax);
    return x;
  }
  double wrapLong(double x) const {
    x += M_PI;               //Move [0,360]
    x = std::fmod(x,2*M_PI); //Map back to [0,360]
    if(x<0)
      x += 2*M_PI;
    x -= M_PI;               //Move back to [-180,180] system
    return x;
  }
  std::uniform_int_distribution<> coord_dist;
  std::normal_distribution<> mut_dist;
  int fail_max   = 20;
  int fail_count = 0;
  int getNextCoord() {
    return coord_dist(rand_engine());
  }
  double getMutation() {
    return mut_dist(rand_engine());
  }
  Point2D mutateBest(){
    Point2D mutated = best.pt;
    int nextcoord = getNextCoord();
    switch(nextcoord){
      case 0:
        mutated.x += getMutation();
        mutated.x  = wrapLong(mutated.x);
        break;
      case 1:
        mutated.y += getMutation();
        mutated.y  = wrapPolar(mutated.y,-M_PI/2,M_PI/2);
        break;
      default:
        throw std::runtime_error("Unrecognized coordinate to mutate!");
    }
    return mutated;
  }
 public:
  HillClimber(Point2D start0, int fail_max0, double mutation_std){
    coord_dist = std::uniform_int_distribution<>(0,coords-1);
    mut_dist   = std::normal_distribution<>(0,mutation_std);
    start      = PointWithStats(0,start0,"");
    best       = start;
    fail_max   = fail_max0;
  }
  void climb(const PointCloud &wgs84pc){
    best.dist = DistanceToCoast(best.pt,wgs84pc);
    while(fail_count<fail_max){
      steps++;
      const Point2D cand_coord = mutateBest();
      const auto    cand       = PointWithStats(DistanceToCoast(cand_coord,wgs84pc), cand_coord, "");
      if(cand.dist<best.dist){
        fail_count++;
        continue;
      }
      best       = cand;
      fail_count = 0;
    }
  }
  void reset(){
    best       = start;
    fail_count = 0;
    steps      = 0;
  }
};



PointWithStats GetBestPoint(const std::vector<PointWithStats> &pts){
  PointWithStats best = pts.front();
  for(auto &x: pts)
    if(x.dist>best.dist)
      best = x;
  return best;
}



//Run `attempts` different hillclimbing runs, each of which gives up after
//`fail_max` attempts at improvement
PointWithStats HillClimb(
  const Point2D          &pt,
  const PointCloud       &wgs84pc,
  const int    attempts,
  const int    fail_max,
  const double mutation_std
){
  PointWithStats pws = PointWithStats(DistanceToCoast(pt,wgs84pc),pt,"");
  std::vector<PointWithStats> bestv(omp_get_max_threads(), pws);

  HillClimber hc(pt,fail_max,mutation_std);

  //Start a large number of hill-climbing walks from the origin
  #pragma omp parallel for default(none) schedule(static) firstprivate(hc) shared(bestv,wgs84pc,attempts)
  for(int i=0;i<attempts;i++){
    hc.reset();
    hc.climb(wgs84pc);
    if(hc.best.dist>bestv.at(omp_get_thread_num()).dist)
      bestv.at(omp_get_thread_num()) = hc.best;
  }

  return GetBestPoint(bestv);
}



PointWithStats ComplexHillClimb(
  const Point2D &pt,
  const PointCloud &wgs84pc
){
  auto best = HillClimb(pt,     wgs84pc,4*50,  20, 0.3*DEG_TO_RAD);
  best      = HillClimb(best.pt,wgs84pc,4*50,  50, 0.1*DEG_TO_RAD);
  best      = HillClimb(best.pt,wgs84pc,24*50,100,0.05*DEG_TO_RAD);
  return best;
}



std::vector<PointWithStats> FilterPoints(const std::vector<PointWithStats> &pts){
  //Build a point cloud
  PointCloud ptcloud;
  for(const auto &p: pts)
    ptcloud.addPoint(TransLLto3D(p.pt));
  ptcloud.buildIndex();

  std::vector<bool> dominated(pts.size(),false);

  for(unsigned int i=0;i<pts.size();i++){
    const auto neighbors = ptcloud.queryByDistance(TransLLto3D(pts[i].pt),20);
    for(const auto &n: neighbors){
      if(pts.at(i).dist>pts.at(n.first).dist)
        dominated.at(n.first) = true;
    }
  }

  std::vector<PointWithStats> undom;
  for(unsigned int i=0;i<pts.size();i++)
    if(!dominated.at(i))
      undom.push_back(pts.at(i));

  return undom;
}


TEST_CASE("Test with data [expensive]"){
  SetupForProjection("ellipsoidal");

  SUBCASE("PointInLandmass"){
    std::vector<Point2D> cities;
    cities.emplace_back(95.0175,27.4833); //MOHANBARI, INDIA
    cities.emplace_back(-73.0603,-40.6114); //OSORNO, CHILE
    cities.emplace_back(92.9786,24.9128); //SILCHAR, INDIA
    cities.emplace_back(22.3931,62.4625); //KAUHAJOKI, FINLAND
    cities.emplace_back(14.4228,62.0478); //SVEG, SWEDEN
    cities.emplace_back(-81.755,26.5361); //FORT MYERS, USA
    cities.emplace_back(150.988,-33.9244); //SYDNEY, AUSTRALIA
    cities.emplace_back(-1.09556,52.0408); //TURWESTON, U.K.
    cities.emplace_back(22.0189,50.11); //RZESZOW, POLAND
    cities.emplace_back(-77.5983,-9.34722); //ANTA, PERU
    cities.emplace_back(-1.1775,47.4081); //ANCENIS, FRANCE
    cities.emplace_back(8.39889,46.9747); //BUOCHS, SWITZERLAND
    cities.emplace_back(-154.911,59.7536); //ILIAMNA, USA
    cities.emplace_back(-77.7956,24.6978); //ANDROS TOWN, BAHAMAS
    cities.emplace_back(-111.348,25.9892); //LORETO, MEXICO
    cities.emplace_back(25.8225,-17.8217); //LIVINGSTONE, ZAMBIA
    cities.emplace_back(-71.3439,-17.695); //ILO, PERU
    cities.emplace_back(-1.52306,43.4683); //BIARRITZ-BAYONNE, FRANCE
    cities.emplace_back(39.0081,8.71556); //DEBRE ZEIT, ETHIOPIA
    cities.emplace_back(-59.2278,-37.2372); //TANDIL, ARGENTINA
    cities.emplace_back(124.611,8.41444); //LADAG, PHILIPPINES
    cities.emplace_back(86.1486,23.6433); //BOKARO, INDIA
    cities.emplace_back(-89.0558,13.4406); //SAN SALVADOR, EL SALVADOR
    cities.emplace_back(30.3986,-29.6489); //PIETERMARITZBURG, SOUTH AFRICA
    cities.emplace_back(-48.4761,-1.37917); //BELEM, BRAZIL
    cities.emplace_back(42.5858,16.9011); //GIZAN, SAUDI ARABIA
    cities.emplace_back(25.6933,61.1439); //VESIVEHMAA, FINLAND
    cities.emplace_back(48.6325,-13.4847); //AMPAMPAMENA, MADAGASCAR
    cities.emplace_back(2.85944,30.5711); //EL GOLEA, ALGERIA
    cities.emplace_back(14.1372,57.2922); //HAGSHULT, SWEDEN
    cities.emplace_back(-50.6514,-24.3175); //TELEMACO BORBA, BRAZIL
    cities.emplace_back(15.4394,46.9908); //GRAZ, AUSTRIA
    cities.emplace_back(94.9292,21.1819); //BAGAN, MYANMAR
    cities.emplace_back(-86.6847,34.6786); //REDSTONE, USA
    cities.emplace_back(3.42361,46.5344); //MOULINS, FRANCE
    cities.emplace_back(-4.96111,51.8331); //HAVERFORDWEST, ENGLAND
    cities.emplace_back(-61.7892,50.1897); //NATASHQUAN, CANADA
    cities.emplace_back(-81.3894,35.7411); //HICKORY, USA
    cities.emplace_back(-3.50333,39.9375); //OCANA, SPAIN
    cities.emplace_back(-84.5872,42.7786); //LANSING, USA
    cities.emplace_back(101.743,6.51972); //NARATHIWAT, THAILAND
    cities.emplace_back(-0.00638889,43.1786); //TARBES, FRANCE
    cities.emplace_back(33.6247,34.875); //LARNACA, CYPRUS
    cities.emplace_back(21.3097,45.1467); //VRSAC, YUGOSLAVIA
    cities.emplace_back(-103.603,35.1828); //TUCUMCARI, USA
    cities.emplace_back(-94.3672,35.3364); //FORT SMITH, USA
    cities.emplace_back(19.3981,51.7219); //LODZ, POLAND
    cities.emplace_back(14.1875,48.2331); //LINZ, AUSTRIA
    cities.emplace_back(128.881,27.8361); //TOKUNOSHIMA, JAPAN
    cities.emplace_back(7.66722,53.5478); //WITTMUNDHAFEN, GERMANY
    for(auto &x: cities)
      x.toRadians();

    for(unsigned int i=0;i<cities.size();i++)
    for(unsigned int j=0;j<cities.size();j++)
      CHECK(std::abs(GeoDistanceEllipsoid(cities.at(i),cities.at(j))-GeoDistanceSphere(cities.at(i),cities.at(j)))<=37);

    for(const auto &c: cities){
      const auto converted = EllipsoidCartesiantoWGS84(WGS84toEllipsoidCartesian(c));
      CHECK(c.x==doctest::Approx(converted.x));
      CHECK(c.y==doctest::Approx(converted.y));
    }

    for(const auto &c: cities){
      const auto converted = SphericalCartesiantoWGS84(WGS84toSphericalCartesian(c));
      CHECK(c.x==doctest::Approx(converted.x));
      CHECK(c.y==doctest::Approx(converted.y));
    }
  }
}

void PrintPoints(
  const std::string &coastname,
  const unsigned int poi_num,
  const PointCloud  &wgs84pc,
  PointWithStats     pt,
  std::ofstream     &fout,
  std::ofstream     &fout_circ
){
  const auto pt3d = TransLLto3D(pt.pt);

  //Get the three nearest points. These define the circle of inaccessibility.
  //Points are returned sorted by distance
  const auto circ_pts = wgs84pc.kNN(pt3d, 3);

  //Get the farthest of the three (should be within epsilon of the nearest)
  const auto dist_to_farthest = std::sqrt(circ_pts.at(2).second);

  //auto annulus = wgs84pc.queryByAnnulus(pt3d, dist_to_farthest, dist_to_farthest+1);
  auto annulus = wgs84pc.queryByDistance(pt3d, dist_to_farthest+30);

  // ClockwiseSorter(annulus, pt3d, wgs84pc);

  if(annulus.size()<3){
    std::cerr<<"Warning: Point "<<pt.pt.x<<" "<<pt.pt.y<<" had no annulus. Skipping."<<std::endl;
    return;
  }

  // std::ostringstream geom;
  // geom << "POLYGON((";

  for(unsigned int i=0;i<annulus.size();i++){
    auto cp = Trans3DtoLL(wgs84pc.pts.at(annulus.at(i).first));
    cp.toDegrees();
    fout_circ<<poi_num<<","
             <<i      <<","
             <<coastname        <<","
             <<chosen_projection<<","
             <<std::fixed<<std::setprecision(10)<<cp.x<<","
             <<std::fixed<<std::setprecision(10)<<cp.y<<","
             <<std::fixed<<std::setprecision(10)<<GeoDistance(pt.pt,cp)
             <<"\n";
    // geom<<std::fixed<<std::setprecision(10)<<cp.x<<" ";
    // geom<<std::fixed<<std::setprecision(10)<<cp.y<<",";
  }
  // {
  //   const auto frontpt = Trans3DtoLL(wgs84pc.pts.at(annulus.at(0).first));
  //   geom<<frontpt.x<<" "<<frontpt.y;
  // }
  // geom<<"))";

  pt.pt.toDegrees();
  fout<<poi_num<<","
      <<coastname<<","
      <<chosen_projection<<","
      <<std::fixed<<std::setprecision(10)<<pt.pt.x<<","
      <<std::fixed<<std::setprecision(10)<<pt.pt.y<<","
      <<std::fixed<<std::setprecision(10)<<pt.dist<<","
      <<pt.label<<","
      //<<"\""<<geom.str()<<"\""
      <<"\"\""
      <<"\n";
}



#ifdef DOCTEST_CONFIG_DISABLE

int main(int argc, char **argv){
  CLI::App app{"Pole of Inaccesibility finder"};

  std::string coastname;
  std::string shapefile_name;
  std::string shapefile_layer_name;
  std::string projection_name;
  std::string out_prefix;
  std::string previous_poles_fname;
  double initial_point_spacing = 300;
  double south_of = 90;
  double north_of = -90;
  double interp_dist = 0.5;
  double distance_filter = -std::numeric_limits<double>::infinity();
  double likely_frac = 1;
  app.add_option("-c,--coastname",    coastname,                 "Coastline name to use in the output CSV")->required();
  app.add_option("-s,--shapefile",    shapefile_name,            "Shapefile to use as input data source")->required();
  app.add_option("-l,--layer",        shapefile_layer_name,      "Layer in the shapefile to use as input data source")->required();
  app.add_option("-p,--projection",   projection_name,           "Which projection to use: spherical, ellipsoidal, haversine")->required();
  app.add_option("-o,--output",       out_prefix,                "Prefix to use for output (Default: poi-<coastname>-<projection>)");
  app.add_option("--spacing",         initial_point_spacing,     "Desired spacing of initial points in kilometers (Default: 300)");
  app.add_option("--south_of",        south_of,                  "Generated points must be south of this latitude (in degrees) (Default: 90)");
  app.add_option("--north_of",        north_of,                  "Generated points must be north of this latitude (in degrees) (Default: -90)");
  app.add_option("--interp_dist",     interp_dist,               "Interpolate edges so that no edge is longer than this in kilometers (Default: 0.5)");
  app.add_option("--distance_filter", distance_filter,           "Filter out initial poles closer to the coastline than this (in km) (Default: -inf)");
  app.add_option("--likely_frac",     likely_frac,               "Proportion of initial points to consider likely for further exploration (in [0,1]) (Default: 1)");
  app.add_option("--previous_poles",  previous_poles_fname,      "Space-delimited file of previous pole locations to try to improve on (Lng Lat Label)");

  CLI11_PARSE(app, argc, argv);

  SetupForProjection(projection_name);

  if(out_prefix.empty()){
    out_prefix = "poi-"+coastname+"-"+projection_name;
  }

  //PREVIOUS POLES
  std::vector<PointWithStats> previous_poles;
  if(!previous_poles_fname.empty()){
    std::ifstream fin(previous_poles_fname);
    double lng;
    double lat;
    std::string label;
    std::cerr<<"Loading previous points (Lng Lat Label)..."<<std::endl;
    while(fin>>lng>>lat){
      std::getline(fin, label);
      previous_poles.emplace_back(99999.0, Point2D{lng, lat}, "prev," + label);
      std::cerr<<"\t"<<lng<<" "<<lat<<" "<<label<<std::endl;
    }
  }

  PointCloud wgs84pc;
  wgs84pc = ReadPointCloudFromShapefile(shapefile_name, shapefile_layer_name, interp_dist);
  wgs84pc.buildIndex();

  std::cerr<<"Getting distance to previous poles..."<<std::endl;
  for(auto &pp: previous_poles){
    pp.pt.toRadians();
    pp.dist = DistanceToCoast(pp.pt, wgs84pc);
  }

  std::ofstream fout(out_prefix + "-poles_of_inaccessibility.csv");
  std::ofstream fout_circ(out_prefix + "-circles_of_inaccessibility.csv");
  fout<<"poi_num,data,proj,PoleX,PoleY,Distance,Type,Label,geom\n";
  fout_circ<<"poi_num,circle_pt_num,data,proj,X,Y,distance\n";

  //Generate evenly-spaced points covering the whole globe
  OrientationGenerator og(initial_point_spacing, south_of, north_of);

  std::vector<PointWithStats> extrema;

  std::cerr<<"Finding likely initial points..."<<std::endl;
  ProgressBar pg(og.size());
  for(int i=0;i<og.size();i++){
    //Get the next initial pole
    const auto initial_pole = og(i);

    //Ignore poles outside of the bounds of interest
    if(!initial_pole.first){
      continue;
    }

    const auto distance = DistanceToCoast(initial_pole.second, wgs84pc);

    // Ignore poles that are distance-filtered
    if(distance <= distance_filter){
      continue;
    }

    extrema.emplace_back(distance, initial_pole.second, "extrema");
    ++pg;
  }
  pg.stop();
  std::cerr<<"Found "<<extrema.size()<<" likely initial points."<<std::endl;

  //Retain only the top `likely_frac` of points
  std::sort(extrema.begin(), extrema.end(), [](const auto &a, const auto &b) { return a.dist > b.dist; });
  const size_t cut_point = static_cast<size_t>(std::ceil(likely_frac * extrema.size()));
  extrema.resize(cut_point);
  extrema.shrink_to_fit();
  std::cerr<<"Considering "<<extrema.size()<<" of the initial points for hill climbing."<<std::endl;
  std::cerr<<"Most isolated point is at "<<extrema.front().dist<<"km and least isolated is at "<<extrema.back().dist<<"km."<<std::endl;

  std::cerr<<"Performing hill climbing from likely initial points..."<<std::endl;
  pg = ProgressBar(extrema.size());
  for(auto &initial_pole: extrema){
    initial_pole = ComplexHillClimb(initial_pole.pt, wgs84pc);
    ++pg;
  }
  pg.stop();

  std::cout<<"Filtering points to reduce density..."<<std::endl;
  extrema = FilterPoints(extrema);

  //Label extrema to indicate origin of points
  for(auto &x: extrema){
    x.label = "ours," + x.label;
  }

  std::cerr<<"Performing hill climbing from previous points..."<<std::endl;
  pg = ProgressBar(20*previous_poles.size());
  for(const auto &pp: previous_poles)
  for(int i=0;i<20;i++){
    extrema.push_back(ComplexHillClimb(pp.pt, wgs84pc));
    extrema.back().label = "improved+" + pp.label;
    ++pg;
  }
  pg.stop();

  std::cerr<<"Sorting all the extrema we found..."<<std::endl;
  std::sort(extrema.begin(), extrema.end(), [](const auto &a, const auto &b) { return a.dist > b.dist; });
  std::cerr<<"Most isolated point is at "<<extrema.front().dist<<"km and least isolated is at "<<extrema.back().dist<<"km."<<std::endl;

  for(size_t i=0;i<extrema.size();i++){
    PrintPoints(coastname, i, wgs84pc, extrema.at(i), fout, fout_circ);
  }

  return 0;
}

#endif