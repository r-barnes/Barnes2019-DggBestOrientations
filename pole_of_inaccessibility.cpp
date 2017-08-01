#ifndef DOCTEST_CONFIG_DISABLE
  #define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#endif

#include "Polygon.hpp"
#include "SpIndex.hpp"
#include "PointCloud.hpp"
#include "Solid.hpp"
#include "GeoStuff.hpp"
#include "Point.hpp"
#include "Orientation.hpp"
#include "Timer.hpp"
#include "IndexedShapefile.hpp"
#include "OrientationIndex.hpp"
#include "Progress.hpp"
#include "random.hpp"
#include <iostream>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <array>
#include <iomanip>
#include <cassert>
#include <fstream>
#include <functional>
#include "doctest.h"
#include <omp.h>

#ifdef ENV_XSEDE
  const std::string FILE_WGS84_LANDMASS = "/home/rbarnes1/scratch/dgg_best/land-polygons-complete-4326/land_polygons.shp";
  const std::string FILE_OUTPUT_PREFIX  = "/home/rbarnes1/scratch/dgg_best/";
  const std::string FILE_MERC_LANDMASS  = "/home/rbarnes1/scratch/dgg_best/land-polygons-split-3857/land_polygons.shp";
#elif ENV_CORI
  const std::string FILE_WGS84_LANDMASS = "/global/homes/r/rbarnes/dgg_best/land-polygons-complete-4326/land_polygons.shp";
  const std::string FILE_OUTPUT_PREFIX  = "/global/homes/r/rbarnes/dgg_best/";
  const std::string FILE_MERC_LANDMASS  = "/global/homes/r/rbarnes/dgg_best/land-polygons-split-3857/land_polygons.shp";
#elif ENV_LAPTOP
  const std::string FILE_WGS84_LANDMASS = "data/land-polygons-complete-4326/land_polygons.shp";
  const std::string FILE_OUTPUT_PREFIX  = "/z/";
  const std::string FILE_MERC_LANDMASS  = "data/land-polygons-split-3857/land_polygons.shp";
#else
  #error "ENV_XSEDE or ENV_LAPTOP must be defined!"
#endif

const double Rearth = 6371; //km

const double DEG_TO_RAD = M_PI/180.0;
const double RAD_TO_DEG = 180.0/M_PI;
const double MAX_COAST_INTERPOINT_DIST = 0.5; //km
const double EDGE_OVERLAPS_SAMPLE_DIST = 10; //km

//Polyhedron and Projection configuration globals
TransLLto3D_t TransLLto3D;
Trans3DtoLL_t Trans3DtoLL;
GeoDistance_t GeoDistance;
std::string chosen_projection;

typedef std::pair<double, Point2D> PointWithStats;


void SetupForProjection(const std::string projection){
  chosen_projection = projection;
  if(projection=="spherical"){
    TransLLto3D = WGS84toSphericalCartesian;
    Trans3DtoLL = SphericalCartesiantoWGS84;
    GeoDistance = GeoDistanceSphere;
  } else if(projection=="ellipsoidal"){
    TransLLto3D = WGS84toEllipsoidCartesian;
    Trans3DtoLL = EllipsoidCartesiantoWGS84;
    GeoDistance = GeoDistanceEllipsoid;
  } else {
    throw std::runtime_error("Unrecognized projection!");
  }
}



PointCloud ReadPointCloudFromShapefile(std::string filename, std::string layer){
  std::cerr<<"Reading point cloud from shapefile..."<<std::endl;
  auto landmass_wgs84 = ReadShapefile(filename, layer);

  std::cerr<<"Converting WGS84 polygons to radians..."<<std::endl;
  for(auto &p: landmass_wgs84)
    p.toRadians();

  const auto poly_fill_point = [&](const Point2D &a, const Point2D &b){
    std::vector<Point2D> ret;
    const auto quickdist = GeoDistanceFlatEarth(a,b);
    if(quickdist>MAX_COAST_INTERPOINT_DIST){
      const auto gcg = GreatArcFactory::make(chosen_projection,a,b,MAX_COAST_INTERPOINT_DIST);
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
  #pragma omp parallel for default(none) schedule(static) shared(landmass_wgs84) reduction(merge: wgs84_ll_flat)
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

  PointCloud wgs84pc;
  std::cerr<<"Adding polygon points to PointCloud..."<<std::endl;
  wgs84pc.pts.swap(wgs84_xyz);

  return wgs84pc;
}



//Returns true if the point falls within a landmass
bool PointInLandmass(const Point2D &ll, const IndexedShapefile &landmass){
  if(ll.y>83.7*DEG_TO_RAD) //The island "83-42" as at 83.7N so anything north of this is on water
    return false;
  if(ll.y<-80*DEG_TO_RAD)  //The southmost extent of water is ~79.5S, so anything south of this is on land
    return true;
  auto xy = WGS84toEPSG3857(ll);
  const auto pid = landmass.sp.queryPoint(xy);
  if(pid==-1)
    return false;
  if(landmass.isRect(pid)) //R-tree rectangle is coincident with the bounding box of a rectangular polygon, so no contains point check is necessary
    return true;
  if(landmass.polys.at(pid).containsPoint(xy))
    return true;
  return false;
}



double DistanceToCoast(const Point2D &pt, const PointCloud &wgs84pc, const IndexedShapefile &landmass){
  const auto cp  = wgs84pc.queryPoint(TransLLto3D(pt)); //Closest point
  const auto llc = Trans3DtoLL(cp);
  auto dist = GeoDistance(llc,pt);
  if(PointInLandmass(pt,landmass))
    dist = -dist;
  return dist;
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
    Point2D mutated = best.second;
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
    start      = std::make_pair(0,start0);
    best       = start;
    fail_max   = fail_max0;
  }
  void climb(
    const PointCloud       &wgs84pc,
    const IndexedShapefile &landmass,
    const std::function<bool(const PointWithStats&, const PointWithStats&)> dom_checker
  ){
    best.first = DistanceToCoast(best.second,wgs84pc,landmass);
    while(fail_count<fail_max){
      steps++;
      const Point2D cand_coord = mutateBest();
      const auto    cand       = std::make_pair(DistanceToCoast(cand_coord,wgs84pc,landmass), cand_coord);
      if(!dom_checker(cand,best)){
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



PointWithStats GetBestPoint(
  const std::vector<PointWithStats> &pts,
  const std::function<bool(const PointWithStats&, const PointWithStats&)> dom_checker
){
  PointWithStats best = pts.front();
  for(auto &x: pts)
    if(dom_checker(x,best))
      best = x;
  return best;
}



//Run `attempts` different hillclimbing runs, each of which gives up after
//`fail_max` attempts at improvement
PointWithStats HillClimb(
  const Point2D          &pt,
  const PointCloud       &wgs84pc,
  const IndexedShapefile &landmass,
  const std::function<bool(const PointWithStats&, const PointWithStats&)> dom_checker,
  const int    attempts,
  const int    fail_max,
  const double mutation_std
){
  PointWithStats pws = std::make_pair(DistanceToCoast(pt,wgs84pc,landmass),pt);
  std::vector<PointWithStats> bestv(omp_get_max_threads(), pws);

  HillClimber hc(pt,fail_max,mutation_std);

  //Start a large number of hill-climbing walks from the origin
  #pragma omp parallel for default(none) schedule(static) firstprivate(hc) shared(bestv,landmass,wgs84pc)
  for(int i=0;i<attempts;i++){
    hc.reset();
    hc.climb(wgs84pc, landmass, dom_checker);
    if(dom_checker(hc.best,bestv.at(omp_get_thread_num())))
      bestv.at(omp_get_thread_num()) = hc.best;
  }

  return GetBestPoint(bestv, dom_checker);
}



PointWithStats ComplexHillClimb(
  const Point2D &pt,
  const PointCloud &wgs84pc,
  const IndexedShapefile &landmass,
  const std::function<bool(const PointWithStats&, const PointWithStats&)> dom_checker
){
  auto best = HillClimb(pt,wgs84pc,landmass,dom_checker,24*100,20,0.3*DEG_TO_RAD);
  best      = HillClimb(best.second,wgs84pc,landmass,dom_checker,24*50,100,0.1*DEG_TO_RAD);
  best      = HillClimb(best.second,wgs84pc,landmass,dom_checker,24*50,100,0.05*DEG_TO_RAD);
  best      = HillClimb(best.second,wgs84pc,landmass,dom_checker,24*50,100,0.01*DEG_TO_RAD);
  return best;
}



std::vector<PointWithStats> FilterPoints(
  const std::vector<PointWithStats> &pts,
  const std::function<bool(const PointWithStats&, const PointWithStats&)> dom_checker
){
  //Build a point cloud
  PointCloud ptcloud;
  for(const auto &p: pts)
    ptcloud.addPoint(TransLLto3D(p.second));
  ptcloud.buildIndex();

  std::vector<bool> dominated(pts.size(),false);

  for(unsigned int i=0;i<pts.size();i++){
    const auto neighbors = ptcloud.queryByDistance(TransLLto3D(pts[i].second),5);
    for(const auto &n: neighbors){
      if(dom_checker(pts.at(i),pts.at(n.first)))
        dominated.at(n.first) = true;
      else
        dominated.at(i) = true;
    }
  }

  std::vector<PointWithStats> undom;
  for(unsigned int i=0;i<pts.size();i++)
    if(!dominated.at(i))
      undom.push_back(pts.at(i));

  return undom;  
}


TEST_CASE("Test with data [expensive]"){
  const auto landmass = IndexedShapefile(FILE_MERC_LANDMASS,"land_polygons");

  SetupForProjection("ellipsoidal");

  SUBCASE("Check that Fuller orientation has no overlaps"){
    SolidXY p = OrientationFullerIcosahedron();

    int ocount = 0;
    for(const auto &v: p.v)
      if(PointInLandmass(v,landmass))
        ocount++;

    std::cerr<<"Fuller count: "<<ocount<<std::endl;
    assert(ocount==0);
  }

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

    //Check that all cities are on land
    for(const auto &x: cities)
      CHECK(PointInLandmass(x, landmass)==true);

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



#ifdef DOCTEST_CONFIG_DISABLE

int main(int argc, char **argv){
  if(argc!=2){
    std::cerr<<argv[0]<<" <Projection>"<<std::endl;
    return -1;
  }

  SetupForProjection(argv[1]);

  const auto landmass = IndexedShapefile(FILE_MERC_LANDMASS,"land_polygons");

  PointCloud wgs84pc;
  wgs84pc = ReadPointCloudFromShapefile(FILE_WGS84_LANDMASS, "land_polygons");
  wgs84pc.buildIndex();

  const auto dom_min = [](const PointWithStats &a, const PointWithStats &b){ return a.first<b.first; };
  const auto dom_max = [](const PointWithStats &a, const PointWithStats &b){ return a.first>b.first; };

  //Generate evenly-spaced points covering the whole globe
  OrientationGenerator og(300,180*DEG_TO_RAD);

  std::ofstream fout(FILE_OUTPUT_PREFIX + "poi.csv");

  std::vector<PointWithStats> min_pts;
  std::vector<PointWithStats> max_pts;

  ProgressBar pg(og.size());
  fout<<"Type,PoleX,PoleY,Distance\n";
  for(int i=0;i<og.size();i++){
    auto pt = og(i);
    min_pts.push_back(ComplexHillClimb(pt, wgs84pc, landmass, dom_min));
    max_pts.push_back(ComplexHillClimb(pt, wgs84pc, landmass, dom_max));
  }

  min_pts = FilterPoints(min_pts, dom_min);
  max_pts = FilterPoints(max_pts, dom_min);

  for(auto p: min_pts){
    p.second.toDegrees();
    fout<<"min,"<<p.second.x<<","<<p.second.y<<","<<p.first<<"\n";
  }
  for(auto p: max_pts){
    p.second.toDegrees();
    fout<<"max,"<<p.second.x<<","<<p.second.y<<","<<p.first<<"\n";
  }

  return 0;
}

#endif