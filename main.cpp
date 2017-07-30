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
  const double MUTATION_WIDTH           = 0.01;
#elif ENV_CORI
  const std::string FILE_WGS84_LANDMASS = "/global/homes/r/rbarnes/dgg_best/land-polygons-complete-4326/land_polygons.shp";
  const std::string FILE_OUTPUT_PREFIX  = "/global/homes/r/rbarnes/dgg_best/";
  const std::string FILE_MERC_LANDMASS  = "/global/homes/r/rbarnes/dgg_best/land-polygons-split-3857/land_polygons.shp";
  const double MUTATION_WIDTH           = 0.01;
#elif ENV_LAPTOP
  const std::string FILE_WGS84_LANDMASS = "data/land-polygons-complete-4326/land_polygons.shp";
  const std::string FILE_OUTPUT_PREFIX  = "/z/";
  const std::string FILE_MERC_LANDMASS  = "data/land-polygons-split-3857/land_polygons.shp";
  const double MUTATION_WIDTH           = 0.05;
#else
  #error "ENV_XSEDE or ENV_LAPTOP must be defined!"
#endif

const double Rearth = 6371; //km

const double DEG_TO_RAD = M_PI/180.0;
const double RAD_TO_DEG = 180.0/M_PI;
const double MAX_COAST_INTERPOINT_DIST = 0.5; //km
const double EDGE_OVERLAPS_SAMPLE_DIST = 10; //km

const double FILTER_COARSE_ORIENTATIONS_WITHIN = 100; //km

#ifdef FINE_RESOLUTION //Used for science
  const double COARSE_SPACING      = 100;  //km - Desired interpoint spacing for finding prospective orienations
  const double COARSE_RADIAL_LIMIT = 90*DEG_TO_RAD;
  const double COARSE_THETA_MIN    = 0;
  const double COARSE_THETA_MAX    = 72*DEG_TO_RAD;
  const double COARSE_THETA_STEP   = 1*DEG_TO_RAD;

#elif COARSE_RESOLUTION //Used for profiling
  const double COARSE_SPACING      = 200;  //km - Desired interpoint spacing for finding prospective orienations
  const double COARSE_RADIAL_LIMIT = 90*DEG_TO_RAD;
  const double COARSE_THETA_MIN    = 0;
  const double COARSE_THETA_MAX    = 72*DEG_TO_RAD;
  const double COARSE_THETA_STEP   = 1.0*DEG_TO_RAD;

#else
  this-is-an-error
#endif

const int    OVERLAPS_TO_BEAT    = 8; //Number of overlaps beyond (and including) which we are interested

typedef std::vector< std::vector<unsigned int> > norientations_t;

//Orientations with no vertices on land are always of interest, as are overlaps
//with many orientations on land. This convenience function is used elsewhere to
//filter by overlaps
bool OverlapOfInterest(unsigned char overlap_count){
  return overlap_count==0 || overlap_count>=OVERLAPS_TO_BEAT; //TODO
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
      const GreatCircleGenerator gcg(a,b,MAX_COAST_INTERPOINT_DIST);
      for(unsigned int i=0;i<gcg.size();i++)
        ret.push_back(gcg(i));
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
  #pragma omp parallel for default(none) schedule(static) shared(wgs84_ll_flat) reduction(merge:wgs84_xyz)
  for(unsigned int i=0;i<wgs84_ll_flat.size();i++)
    wgs84_xyz.push_back(WGS84toEllipsoidCartesian(wgs84_ll_flat[i]));
    //wgs84_xyz.push_back(wgs84_ll_flat[i].toXYZ(1));

  wgs84_ll_flat.clear();
  wgs84_ll_flat.shrink_to_fit();

  wgs84_xyz.shrink_to_fit();

  PointCloud wgs84pc;
  std::cerr<<"Adding polygon points to PointCloud..."<<std::endl;
  wgs84pc.pts.swap(wgs84_xyz);

  return wgs84pc;
}



//Determine which orientations are in the local neighbourhood of another
//orientation
template<class T>
norientations_t FindNearbyOrientations(const std::vector<T> &osc, const double distance){
  Timer tmr_bi;
  std::cerr<<"Building kd-tree"<<std::endl;
  OrientationIndex oidx(osc);
  std::cerr<<"Time = "<<tmr_bi.elapsed()<<std::endl;

  std::cerr<<"Finding nearby orientations..."<<std::endl;

  ProgressBar pg(osc.size());
  norientations_t oneighbors(osc.size());
  #pragma omp parallel for default(none) schedule(static) shared(osc,oneighbors,oidx,pg)
  for(unsigned int i=0;i<osc.size();i++){
    oneighbors[i] = oidx.query(i,distance);
    ++pg;
  }

  std::cerr<<"Time taken = "<<pg.stop()<<std::endl;

  return oneighbors;
}



OCollection FilterOutOrienationsWithNeighbours(const OCollection &oc, const double distance){
  const auto norients = FindNearbyOrientations(oc, distance);
  std::vector<bool> keep(oc.size(),true);
  for(unsigned int focal=0;focal<norients.size();focal++){
    if(!keep.at(focal))
      continue;
    for(const auto &n: norients[focal])
      keep.at(n) = false;
  }

  OCollection keep_these;
  for(unsigned int o=0;o<oc.size();o++)
    if(keep.at(o))
      keep_these.push_back(oc.at(o));

  return keep_these;
}

OSCollection FilterOutDominatedOrienations(
  const OSCollection &osc,
  const double distance,
  std::function<bool(const OrientationWithStats&, const OrientationWithStats&)> dom_checker
){
  std::vector<bool> dominated(osc.size(),false);
  const auto norients = FindNearbyOrientations(osc, distance);

  for(unsigned int focal=0;focal<norients.size();focal++){
    for(const auto &n: norients[focal])
      if(dom_checker(osc[focal],osc[n])) //If focal dominats n
        dominated.at(n) = true;
  }

  OSCollection keep_these;
  for(unsigned int o=0;o<osc.size();o++)
    if(!dominated.at(o))
      keep_these.push_back(osc.at(o));

  return keep_these;
}



std::ostream& operator<<(std::ostream &out, const std::vector<bool> &v){
  for(const auto x: v)
    out<<(int)x;
  return out;
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



//Returns a bitset indicating indicating which vertices of a polyhedron lie
//within a landmass
auto OrientationOverlaps(const SolidXY &sxy, const IndexedShapefile &landmass){
  std::vector<bool> overlaps(sxy.v.size(),false);
  for(unsigned int vi=0;vi<sxy.v.size();vi++)
    if(PointInLandmass(sxy.v[vi],landmass))
      overlaps.at(vi) = true;
  return overlaps;
}



//Given a set of orientations, return a new set in which only those orientations
//which fall entirely in the water or have a certain number of vertices on land
//are returned in a new set. Orientations are not returned in the same order.
OCollection OrientationsFilteredByOverlaps(
  const IndexedShapefile &landmass
){
  std::cerr<<"Filtering orientations for overlaps..."<<std::endl;

  const OrientationGenerator og(COARSE_SPACING, COARSE_RADIAL_LIMIT);

  std::cerr<<"Orientations to generate sans theta-rotation = "<<(long)og.size()<<std::endl;

  const long orienations_to_search = ((long)(og.size()*((COARSE_THETA_MAX-COARSE_THETA_MIN)/COARSE_THETA_STEP)));
  std::cerr<<"Orientations to generate with theta-rotation = "<<orienations_to_search<<std::endl;

  ProgressBar pg(og.size());
  OCollection ret;
  #pragma omp declare reduction (merge : OCollection : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
  #pragma omp parallel for default(none) schedule(static) shared(landmass,pg) reduction(merge: ret)
  for(unsigned int o=0;o<og.size();o++){
    const auto pole = og(o);
    for(double theta=COARSE_THETA_MIN;theta<=COARSE_THETA_MAX;theta+=COARSE_THETA_STEP){
      const Orientation ori(pole,theta);
      SolidXY sxy(ori);
      const auto overlaps     = OrientationOverlaps(sxy, landmass);
      const int overlap_count = std::accumulate(overlaps.begin(),overlaps.end(),0);
      if(OverlapOfInterest(overlap_count))
        ret.push_back(ori);
    }
    ++pg;
  }

  std::cout<<"Filtering: Time taken = " <<pg.stop()  <<"s"<<std::endl;
  std::cerr<<"Filtering: Found = "      <<ret.size() <<std::endl;

  return ret;
}



//For a great circle connecting two points, generate sample points along the
//circle. Then count how many of these sample points fall within landmasses.
//Maximum value returned is `num_pts+1`
unsigned int GreatCircleOverlaps(const IndexedShapefile &landmass, const Point2D &a, const Point2D &b, const double spacing){
  GreatCircleGenerator gcg(a,b,spacing);

  unsigned int edge_overlaps = 0;
  for(unsigned int i=0;i<=gcg.size();i++)
    edge_overlaps += PointInLandmass(gcg(i),landmass);

  return edge_overlaps;
}



//For all the great circle edges of a polyhedron, determine how many sample
//points along the circles fall within landmasses.
unsigned int OrientationEdgeOverlaps(const SolidXY &sxy, const IndexedShapefile &landmass){
  static const auto   neighbors = SolidXY().neighbors(); //Get a list of neighbouring vertices on the polyhedron
  unsigned int edge_overlaps = 0;
  for(unsigned int n=0;n<neighbors.size();n+=2){
    const auto &a = sxy.v[neighbors[n]];
    const auto &b = sxy.v[neighbors[n+1]];
    edge_overlaps += GreatCircleOverlaps(landmass, a, b, EDGE_OVERLAPS_SAMPLE_DIST);
  }
  return edge_overlaps;
}



//Generate distance statistic for an orientation
OrientationWithStats OrientationStats(const Orientation &o, const PointCloud &wgs84pc, const IndexedShapefile &landmass, const bool do_edge){
  OrientationWithStats ows(o);
  SolidXY sxy(o);

  ows.overlaps      = OrientationOverlaps(sxy, landmass);
  ows.overlap_count = std::accumulate(ows.overlaps.begin(),ows.overlaps.end(),0);

  for(unsigned int i=0;i<sxy.v.size();i++){
    const auto cp  = wgs84pc.queryPoint(WGS84toEllipsoidCartesian(sxy.v[i])); //Closest point
    const auto llc = EllipsoidCartesiantoWGS84(cp);
    auto dist      = GeoDistanceEllipsoid(llc,sxy.v[i]);
    //auto dist    = GeoDistanceFlatEarth(llc,sxy.v[i]);
    if(ows.overlaps.at(i))
      dist = -dist;
    ows.mindist = std::min(ows.mindist,dist);
    ows.maxdist = std::max(ows.maxdist,dist);
    ows.avgdist += dist;
  }

  if(do_edge)
    ows.edge_overlaps = OrientationEdgeOverlaps(sxy, landmass);

  return ows;
}



class HillClimber {
 public:
  static const int coords = 3;
  int steps        = 0;
  OrientationWithStats start;
  OrientationWithStats best;
 private:
  double wrapPolar(double x, const double rangemin, const double rangemax) const {
    if(x<rangemin)
      x = std::abs(x-rangemin)+rangemin;
    else if(x>rangemax)
      x = rangemax-std::abs(x-rangemax);
    return x;    
  }
  double wrapLong(double x, const double rangemin, const double rangemax) const {
    if(x<rangemin)
      x = rangemax + (x-rangemin);
    else if(x>rangemax)
      x = rangemin + (x-rangemax);
    return x;
  }
  std::uniform_int_distribution<> coord_dist;
  std::normal_distribution<> mut_dist;
  int fail_max = 20;
  int fail_count = 0;
  int getNextCoord() {
    return coord_dist(rand_engine());
  }
  double getMutation() {
    return mut_dist(rand_engine());
  }
  Orientation mutateBest(){
    Orientation mutated = best;
    int nextcoord = getNextCoord();
    switch(nextcoord){
      case 0:
        mutated.pole.x += getMutation();
        mutated.pole.x = wrapLong(mutated.pole.x,-M_PI,M_PI);
        break;
      case 1:
        mutated.pole.y += getMutation();
        mutated.pole.y = wrapPolar(mutated.pole.y,-M_PI/2,M_PI/2);
        break;
      case 2:
        mutated.theta += getMutation();
        mutated.theta = wrapLong(mutated.theta,-M_PI,M_PI);
        break;
      default:
        throw std::runtime_error("Unrecognized coordinate to mutate!");
    }
    return mutated;
  }
 public:
  HillClimber(OrientationWithStats start0, int fail_max0, double mutation_std){
    coord_dist = std::uniform_int_distribution<>(0,coords-1);
    mut_dist   = std::normal_distribution<>(0,mutation_std);
    start      = start0;
    best       = start;
    fail_max   = fail_max0;
  }
  void climb(
    const PointCloud &wgs84pc,
    const IndexedShapefile &landmass,
    bool do_edge,
    std::function<bool(const OrientationWithStats&, const OrientationWithStats&)> dom_checker
  ){
    while(fail_count<fail_max){
      steps++;
      auto cand_orient    = mutateBest();
      auto cand_orient_ws = OrientationStats(cand_orient, wgs84pc, landmass, do_edge);
      if(!OverlapOfInterest(cand_orient_ws.overlap_count)){
        fail_count++;
        continue;
      }
      if(!dom_checker(cand_orient_ws,best)){
        fail_count++;
        continue;
      }
      best       = cand_orient_ws;
      fail_count = 0;
    }
  }
  void reset(){
    best       = start;
    fail_count = 0;
    steps      = 0;
  }
};



OrientationWithStats GetBestOrientationWithStats(
  const OSCollection &osc,
  std::function<bool(const OrientationWithStats&, const OrientationWithStats&)> dom_checker
){
  OrientationWithStats best = osc.front();
  for(auto &x: osc)
    if(dom_checker(x,best))
      best = x;
  return best;
}



//Run `attempts` different hillclimbing runs, each of which gives up after
//`fail_max` attempts at improvement
OrientationWithStats HillClimb(
  const Orientation &orient,
  const PointCloud &wgs84pc,
  const IndexedShapefile &landmass,
  bool do_edge,
  std::function<bool(const OrientationWithStats&, const OrientationWithStats&)> dom_checker,
  const int attempts,
  const int fail_max,
  const double mutation_std
){
  OrientationWithStats ows = OrientationStats(orient, wgs84pc, landmass, do_edge);
  std::vector<OrientationWithStats> bestv(omp_get_max_threads(), ows);

  HillClimber hc(ows,fail_max,mutation_std);

  //Start a large number of hill-climbing walks from the origin
  #pragma omp parallel for default(none) schedule(static) firstprivate(hc) shared(bestv,landmass,wgs84pc,do_edge,dom_checker)
  for(int i=0;i<attempts;i++){
    hc.reset();
    hc.climb(wgs84pc, landmass, do_edge, dom_checker);
    if(dom_checker(hc.best,bestv.at(omp_get_thread_num())))
      bestv.at(omp_get_thread_num()) = hc.best;
  }

  return GetBestOrientationWithStats(bestv, dom_checker);
}



OrientationWithStats ComplexHillClimb(
  const Orientation &orient,
  const PointCloud &wgs84pc,
  const IndexedShapefile &landmass,
  bool do_edge,
  std::function<bool(const OrientationWithStats&, const OrientationWithStats&)> dom_checker
){
  auto best = HillClimb(orient,wgs84pc,landmass,do_edge,dom_checker,1000,20,0.1*DEG_TO_RAD);
  best      = HillClimb(best,wgs84pc,landmass,do_edge,dom_checker,48,100,0.1*DEG_TO_RAD);
  return best;
}



std::ofstream& PrintPOI(std::ofstream& fout, const OSCollection &osc, const int i, bool header){
  if(header){
    fout<<"num"           <<","
        <<"overlaps"      <<","
        <<"overlaps"      <<","
        <<"pole.y"        <<","
        <<"pole.x"        <<","
        <<"rtheta"        <<","
        <<"mindist"       <<","
        <<"maxdist"       <<","
        <<"avgdist"       <<","
        <<"edge_overlaps"
        <<"\n";
  } else {
    fout<<i<<","
        <<osc[i].overlaps            <<","
        <<osc[i].overlap_count       <<","
        <<(osc[i].pole.y*RAD_TO_DEG) <<","
        <<(osc[i].pole.x*RAD_TO_DEG) <<","
        <<(osc[i].theta*RAD_TO_DEG)  <<","
        <<osc[i].mindist             <<","
        <<osc[i].maxdist             <<","
        <<(osc[i].avgdist/12)        <<","
        <<osc[i].edge_overlaps
        <<"\n";
  }
  return fout;
}



std::ofstream& PrintPOICoordinates(std::ofstream& fout, const OSCollection &osc, const int pn, bool header){
  if(header){
    fout<<"Num,Lat,Lon,OnLand\n";
  } else {
    SolidXY p(osc[pn]);
    for(unsigned int i=0;i<p.v.size();i++)
      fout<<pn                    <<","
          <<p.v[i].y*RAD_TO_DEG<<","
          <<p.v[i].x*RAD_TO_DEG<<","
          <<(int)osc[pn].overlaps.at(i)
          <<"\n";
  }
  return fout;
}



void PrintOrientations(
  std::string fileprefix,
  const OSCollection &osc
){
  {
    std::ofstream fout(FILE_OUTPUT_PREFIX + fileprefix + "-rot.csv");
    PrintPOI(fout, osc, 0, true);
    for(unsigned int o=0;o<osc.size();o++)
      PrintPOI(fout, osc, o, false);
  }
  {
    std::ofstream fout(FILE_OUTPUT_PREFIX + fileprefix + "-vert.csv");
    PrintPOICoordinates(fout, osc, 0, true);
    for(unsigned int o=0;o<osc.size();o++)
      PrintPOICoordinates(fout, osc, o, false);
  }
}



OSCollection FindBestHelper(
  const std::string prefix,
  const OCollection &orients,
  const PointCloud &wgs84pc,
  const IndexedShapefile &landmass,
  bool do_edge,
  std::function<bool(const OrientationWithStats&, const OrientationWithStats&)> dom_checker
){
  Timer tmr;
  std::cerr<<"Finding best for '"<<prefix<<"'"<<std::endl;
  OSCollection ret;
  for(unsigned int o=0;o<orients.size();o++){
    std::cerr<<"\t"<<prefix<<" progress="<<o<<"/"<<orients.size()<<", time="<<tmr.elapsed()<<"s, total_est="<<((tmr.elapsed()/o)*orients.size())<<std::endl;
    ret.push_back(ComplexHillClimb(orients.at(o),wgs84pc,landmass,do_edge,dom_checker));
  }

  ret = FilterOutDominatedOrienations(ret,100,dom_checker);

  PrintOrientations(prefix,ret);

  std::cerr<<"Finished '"<<prefix<<"' in time = "<<tmr.elapsed()<<std::endl;

  return ret;
}



void FindBest(
  const OCollection      &orients,
  const PointCloud       &wgs84pc,
  const IndexedShapefile &landmass
){
  Timer tmr;
  std::cerr<<"Find best..."<<std::endl;

  {
    FindBestHelper("out_min_mindist", orients, wgs84pc, landmass, false,
      [](const OrientationWithStats &a, const OrientationWithStats &b){ return a.mindist<b.mindist; }
    );
    FindBestHelper("out_max_mindist", orients, wgs84pc, landmass, false,
      [](const OrientationWithStats &a, const OrientationWithStats &b){ return a.mindist>b.mindist; }
    );
    

    FindBestHelper("out_min_maxdist", orients, wgs84pc, landmass, false,
      [](const OrientationWithStats &a, const OrientationWithStats &b){ return a.maxdist<b.maxdist; }
    );
    FindBestHelper("out_max_maxdist", orients, wgs84pc, landmass, false,
      [](const OrientationWithStats &a, const OrientationWithStats &b){ return a.maxdist>b.maxdist; }
    );

    
    FindBestHelper("out_min_avgdist", orients, wgs84pc, landmass, false,
      [](const OrientationWithStats &a, const OrientationWithStats &b){ return a.avgdist<b.avgdist; }
    );
    FindBestHelper("out_max_avgdist", orients, wgs84pc, landmass, false,
      [](const OrientationWithStats &a, const OrientationWithStats &b){ return a.avgdist>b.avgdist; }
    );
    

    FindBestHelper("out_min_edge_overlaps", orients, wgs84pc, landmass, true,
      [](const OrientationWithStats &a, const OrientationWithStats &b){ return a.edge_overlaps<b.edge_overlaps; }
    );
    FindBestHelper("out_max_edge_overlaps", orients, wgs84pc, landmass, true,
      [](const OrientationWithStats &a, const OrientationWithStats &b){ return a.edge_overlaps>b.edge_overlaps; }
    );
  }
  
  std::cerr<<"Found best in = "<<tmr.elapsed()<<std::endl;
}



TEST_CASE("Check orientation of generated points"){
  //Use NaN to ensure that all points are generated
  OrientationGenerator og(200,180*DEG_TO_RAD);

  std::vector<Point2D> orients;
  for(long i=0;i<og.size();i++)
    orients.push_back(og(i));

  CHECK(orients.front().y>0);
  CHECK(orients.back().y<0);
  std::ofstream fout("test/test_orientations_spiral.csv");
  fout<<"num,lat,lon\n";
  for(unsigned int o=0;o<orients.size();o++)
    fout<<o<<","<<(orients[o].y*RAD_TO_DEG)<<","<<(orients[o].x*RAD_TO_DEG)<<"\n";
}

//Determine the number of orientations in one quadrant of the 3-space
TEST_CASE("Counting orientations [expensive]"){
  const OrientationGenerator og(200,90*DEG_TO_RAD);

  CHECK(og.size()>0);

  int mincount = std::numeric_limits<int>::max();
  int maxcount = std::numeric_limits<int>::lowest();
  Timer tmr;
  #pragma omp parallel for default(none) schedule(static) reduction(min:mincount) reduction(max:maxcount)
  for(long i=0;i<og.size();i++){
    SolidXYZ p = SolidXY(Orientation(og(i),0)).toXYZ(1);
    int count  = 0;
    for(const auto &v: p.v)
      if(v.y>=0 && v.z>=0)
        count++;
    mincount = std::min(count,mincount);
    maxcount = std::max(count,maxcount);
  }
  CHECK(mincount==2);
  CHECK(maxcount==4);
}



TEST_CASE("Test with data [expensive]"){
  const auto landmass = IndexedShapefile(FILE_MERC_LANDMASS,"land_polygons");

  SUBCASE("Check that Fuller orientation has no overlaps"){
    SolidXY p;
    //Fuller's orientation
    p.v = {{
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
    //p.lon = {{10.53620,  -5.24539,  58.15771, 122.3    ,-143.47849, -67.13233, -57.7    ,  36.5215 , 112.86767, 174.7546 ,-121.84229,-169.4638}};
    //p.lat = {{64.7,   2.300882,  10.447378,  39.1,  50.103201,  23.717925, -39.1, -50.1032, -23.717925,  -2.3009, -10.447345, -64.7}};
    p.toRadians();
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
      CHECK(std::abs(GeoDistanceEllipsoid(cities[i],cities[j])-GeoDistanceHaversine(cities[i],cities[j]))<=37);

    for(const auto &c: cities){
      const auto converted = EllipsoidCartesiantoWGS84(WGS84toEllipsoidCartesian(c));
      CHECK(c.x==doctest::Approx(converted.x));
      CHECK(c.y==doctest::Approx(converted.y));
    }
    

    //Check that all GC arcs between cities include at least some land
    for(unsigned int i=0;i<cities.size();i++)
    for(unsigned int j=0;j<cities.size();j++){
      if(i==j)
        continue;
      CHECK(GreatCircleOverlaps(landmass,cities[i],cities[j],100)>0);
    }

    //Route from Minneapolis to Denver should be entirely on land
    {
      auto a = Point2D(-93,45).toRadians();
      auto b = Point2D(-104.9903, 39.7392).toRadians();
      CHECK(GreatCircleOverlaps(landmass,a,b,50)==24); //TODO: Unsure how many points there will be now that we are using distance instead of point counts
    }
  }
}

TEST_CASE("POIindex"){
  OCollection oc;
  oc.emplace_back(Point2D(-93,45).toRadians(), 0);
  oc.emplace_back(Point2D(-93.1,45.1).toRadians(), 0*DEG_TO_RAD);
  oc.emplace_back(Point2D(-93.1,45.1).toRadians(), 72*DEG_TO_RAD);
  oc.emplace_back(Point2D(-93.2,45.1).toRadians(), 0*DEG_TO_RAD);
  oc.emplace_back(Point2D(-93.1,45.2).toRadians(), 0*DEG_TO_RAD);
  oc.emplace_back(Point2D(-93.2,45.2).toRadians(), 0*DEG_TO_RAD);
  oc.emplace_back(Point2D(-93.2,45.2).toRadians(), 36*DEG_TO_RAD);
  oc.emplace_back(Point2D(23,-23.2).toRadians(), 36*DEG_TO_RAD);
  auto oneighbors = FindNearbyOrientations(oc,100);
  CHECK(oneighbors[0].size()==5);
  CHECK(oneighbors[7].size()==0);
}

TEST_CASE("Generate great cicles between points"){
  SolidXY sxy;
  static const auto neighbors = SolidXY().neighbors();             //Get a list of neighbouring vertices on the polyhedron

  std::ofstream fout("test/test_gc_between_verts");
  fout<<"lat,lon\n";

  for(unsigned int n=0;n<neighbors.size();n+=2){
    const auto &a = sxy.v[neighbors[n]];
    const auto &b = sxy.v[neighbors[n+1]];
    const GreatCircleGenerator gcg(a,b,100);
    CHECK(gcg.getSpacing()==100);
    for(unsigned int i=0;i<gcg.size();i++){
      auto temp = gcg(i);
      fout<<(temp.y*RAD_TO_DEG)<<","<<(temp.x*RAD_TO_DEG)<<"\n";
    }
  }
}



#ifdef DOCTEST_CONFIG_DISABLE

int main(){
  std::cout<<"Rearth              = " << Rearth              <<std::endl;
  std::cout<<"COARSE_SPACING      = " << COARSE_SPACING      <<std::endl;
  std::cout<<"COARSE_RADIAL_LIMIT = " << COARSE_RADIAL_LIMIT <<std::endl;
  std::cout<<"COARSE_THETA_MIN    = " << COARSE_THETA_MIN    <<std::endl;
  std::cout<<"COARSE_THETA_MAX    = " << COARSE_THETA_MAX    <<std::endl;
  std::cout<<"COARSE_THETA_STEP   = " << COARSE_THETA_STEP   <<std::endl;


  assert(!omp_get_nested());

  const auto landmass = IndexedShapefile(FILE_MERC_LANDMASS,"land_polygons");

  OCollection orients;
  orients = OrientationsFilteredByOverlaps(landmass);
  std::cout<<"Orientations generated with overlaps of interest = "<<orients.size()<<std::endl;
  orients = FilterOutOrienationsWithNeighbours(orients, FILTER_COARSE_ORIENTATIONS_WITHIN);
  std::cout<<"Orientations remaining after filtering those with neighbours = "<<orients.size()<<std::endl;

  PointCloud wgs84pc;
  wgs84pc = ReadPointCloudFromShapefile(FILE_WGS84_LANDMASS, "land_polygons");
  wgs84pc.buildIndex();

  std::cerr<<"wgs84pc size = "<<wgs84pc.pts.size()<<std::endl;

  //Seed a PRNG for each thread from the machine's entropy
  seed_rand(0);

  std::cerr<<"Climbing hills..."<<std::endl;

  FindBest(orients, wgs84pc, landmass);

  return 0;
}

#endif