//https://gis.stackexchange.com/questions/70800/how-to-read-vector-datasets-using-gdal-library
//https://math.stackexchange.com/questions/44080/approximating-geodesic-distances-on-the-sphere-by-euclidean-distances-of-a-trans
#include "Polygon.hpp"
#include "SpIndex.hpp"
#include "PointCloud.hpp"
#include "Icosa.hpp"
#include "GeoStuff.hpp"
#include "Point.hpp"
#include "POI.hpp"
#include "Timer.hpp"
#include "IndexedShapefile.hpp"
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicLine.hpp>
#include <GeographicLib/Constants.hpp>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <array>
#include <iomanip>
#include <cassert>
#include <fstream>
#include "doctest.h"

#ifdef ENV_XSEDE
  const std::string FILE_WGS84_LANDMASS = "/home/rbarnes1/scratch/dgg_best/land-polygons-complete-4326/land_polygons.shp";
  const std::string FILE_OUTPUT_ROT     = "/home/rbarnes1/scratch/dgg_best/out-rot.csv";
  const std::string FILE_OUTPUT_VERT    = "/home/rbarnes1/scratch/dgg_best/out-vert.csv";
  const std::string FILE_MERC_LANDMASS  = "/home/rbarnes1/scratch/dgg_best/land-polygons-split-3857/land_polygons.shp";
#elif ENV_LAPTOP
  const std::string FILE_WGS84_LANDMASS = "data/land-polygons-complete-4326/land_polygons.shp";
  const std::string FILE_OUTPUT_ROT     = "/z/out-rot.csv";
  const std::string FILE_OUTPUT_VERT    = "/z/out-vert.csv";
  const std::string FILE_MERC_LANDMASS  = "data/land-polygons-split-3857/land_polygons.shp";
#else
  this-is-an-error
#endif

const double Rearth = 6371; //km
const double pspace = 10;  //km - Desired interpoint spacing

const double DEG_TO_RAD = M_PI/180.0;
const double RAD_TO_DEG = 180.0/M_PI;

//1/10th of a degree grid spacing
const double PRECISION  = 0.1;
const double DIV        = 10.0;

//Neighbouring vertices of icosahedron used for generating rotations
const int NA = 0;
const int NB = 2;
const int NC = 4;

//1 degree grid spacing
//const double PRECISION  = 1;  
//const double DIV        = 1;

bool PointOverlaps(const Point2D &ll, const IndexedShapefile &landmass){
  if(ll.y>83.7*DEG_TO_RAD) //The island "83-42" as at 83.7N so anything north of this is on water
    return false;
  if(ll.y<-80*DEG_TO_RAD)  //The southmost extent of water is ~79.5S, so anything south of this is on land
    return true;
  auto xy = WGS84toEPSG3857(ll);
  const auto pid = landmass.sp.queryPoint(xy);
  if(pid==-1)
    return false;
  if(landmass.polys.at(pid).containsPoint(xy))
    return true;
  return false;
}

TEST_CASE("Test with data [data]"){
  const auto landmass = IndexedShapefile(FILE_MERC_LANDMASS,"land_polygons");

  SUBCASE("Check that Fuller orientation has no overlaps"){
    IcosaXY p;
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
      if(PointOverlaps(v,landmass))
        ocount++;

    std::cerr<<"Fuller count: "<<ocount<<std::endl;
    assert(ocount==0);
  }
}



std::vector<Point2D> GenerateOrientations(){
  std::cerr<<"Generating orientations..."<<std::endl;
  std::vector<Point2D> orientations;

  //Number of points to sample
  const int N = (int)(8*M_PI*Rearth*Rearth/std::sqrt(3)/pspace/pspace);

  //Generate orientations
  #pragma omp declare reduction (merge : std::vector<Point2D> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
  #pragma omp parallel for default(none) schedule(static) reduction(merge: orientations)
  for(int i=0;i<N;i++){
    Point2D temp (
      M_PI*(3-std::sqrt(5))*i,
      std::acos(1-(2.0*i+1.0)/N)
    );
    temp.x = std::fmod(temp.x,2*M_PI)-M_PI;
    temp.y = temp.y-M_PI/2;

    if(!(0 <= temp.x && temp.x<=78*DEG_TO_RAD))
      continue;
    if(temp.y<23*DEG_TO_RAD)
      continue;

    orientations.push_back(temp);
  }

  return orientations;
}

//Determine the number of orientations in one quadrant of the 3-space
TEST_CASE("Counting orientations [expensive]"){
  const auto orientations = GenerateOrientations();

  int mincount = std::numeric_limits<int>::max();
  int maxcount = std::numeric_limits<int>::lowest();
  Timer tmr;
  #pragma omp parallel for default(none) schedule(static) reduction(min:mincount) reduction(max:maxcount)
  for(unsigned int oi=0;oi<orientations.size();oi++)
  for(double rtheta=0;rtheta<72.01*DEG_TO_RAD;rtheta+=PRECISION*DEG_TO_RAD){
    IcosaXYZ p = IcosaXY(orientations[oi],rtheta).toXYZ(1);
    int count  = 0;
    for(const auto &v: p.v)
      if(v.y>=0 && v.z>=0)
        count++;
    mincount = std::min(count,mincount);
    maxcount = std::max(count,maxcount);
  }
  CHECK(mincount==2);
  CHECK(maxcount==3);
}



POICollection FindOrientationsOfInterest(const IndexedShapefile &landmass){
  std::cerr<<"Finding poles..."<<std::endl;
  POICollection poic;
  
  const auto orientations = GenerateOrientations();
  
  long count = 0;

  Timer tmr;
  #pragma omp parallel for default(none) schedule(static) shared(poic,std::cerr,landmass) reduction(+:count)
  for(unsigned int oi=0;oi<orientations.size();oi++)
  for(double rtheta=0;rtheta<72.01*DEG_TO_RAD;rtheta+=PRECISION*DEG_TO_RAD){
    count++;
    IcosaXY p(orientations[oi],rtheta);

    std::bitset<12> overlaps = 0;
    for(unsigned int pi=0;pi<p.v.size();pi++)
      if(PointOverlaps(p.v[pi],landmass))
        overlaps.set(pi);
    if(overlaps==0 || overlaps.count()>=8){
      #pragma omp critical
      poic.emplace_back(overlaps,orientations[oi],rtheta);
    }
  }

  double t = tmr.elapsed();
  std::cout << "Time taken = " << t <<"s"<< std::endl;

  std::cerr<<"Checked = "<<count<<std::endl;
  std::cerr<<"Found "<<poic.size()<<" poles of interest."<<std::endl;

  return poic;
}

void DistancesToIcosaXYs(POICollection &poic){
  std::cerr<<"Reading WGS84 shapefile..."<<std::endl;
  Polygons landmass_wgs84 = ReadShapefile(FILE_WGS84_LANDMASS, "land_polygons");

  PointCloud pc;

  for(auto &ll: landmass_wgs84)
    ll.toRadians();

  std::cerr<<"Adding polygon points to PointCloud..."<<std::endl;
  for(const auto &poly: landmass_wgs84)
  for(const auto &ll: poly.exterior)
    pc.addPoint(ll.toXYZ(1));

  std::cerr<<"Building kd-tree..."<<std::endl;
  pc.buildIndex();

  Timer tmr;

  std::cerr<<"Calculating distances to poles..."<<std::endl;
  #pragma omp parallel for default(none) shared(poic,pc)
  for(unsigned int pn=0;pn<poic.size();pn++){
    IcosaXY p(poic[pn].pole,poic[pn].rtheta);

    for(unsigned int i=0;i<p.v.size();i++){
      const auto cp = pc.queryPoint(p.v[i].toXYZ(1)); //Closest point
      auto llc      = cp.toLatLon();
      auto dist     = GeoDistanceFlatEarth(llc,p.v[i]);
      if(poic[pn].overlaps.test(i))
        dist = -dist;
      poic[pn].mindist = std::min(poic[pn].mindist,dist);
      poic[pn].maxdist = std::max(poic[pn].maxdist,dist);
      poic[pn].avgdist += dist;
    }
  }

  std::cout << "Time taken = " << tmr.elapsed() <<"s"<< std::endl;
}

void EdgeOverlaps(const IndexedShapefile &landmass, POICollection &poic){
  Timer tmr;
  std::cerr<<"Calculating edge overlaps..."<<std::endl;
  const GeographicLib::Geodesic geod(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f());
  const auto neighbors = IcosaXY().neighbors();
  const double ndist   = IcosaXY().neighborDistance()*1000;  //Approximate spacing between vertices in metres
  const double spacing = 10e3;                            //Spacing between points = 10km
  const int    num_pts = int(std::ceil(ndist / spacing)); //The number of intervals
  std::cerr<<"Using "<<num_pts<<" with a "<<spacing<<"m spacing to cover "<<ndist<<"m inter-neighbour distance."<<std::endl;
  #pragma omp parallel for default(none) shared(poic,landmass)
  for(unsigned int pn=0;pn<poic.size();pn++){
    IcosaXY p(poic[pn].pole, poic[pn].rtheta);
    for(unsigned int n=0;n<neighbors.size();n+=2){
      const auto &a = p.v[neighbors[n]];
      const auto &b = p.v[neighbors[n+1]];
      const GeographicLib::GeodesicLine line = geod.InverseLine(
        a.y*RAD_TO_DEG,
        a.x*RAD_TO_DEG,
        b.y*RAD_TO_DEG,
        b.x*RAD_TO_DEG
      );
      const double da = line.Arc() / num_pts;
      for(int i=0;i<=num_pts;i++) {
        Point2D temp;
        line.ArcPosition(i * da, temp.y, temp.x);
        temp.toRadians();
        poic[pn].edge_overlaps += PointOverlaps(temp,landmass);
      }
    }
  }
  std::cout << "Time taken = " << tmr.elapsed() <<"s"<< std::endl;
}

//CHEESE
template<class T>
std::vector<size_t> Dominants(
  const POICollection &poic,
  T dom_checker
){
  std::vector<size_t> dominates(poic.size());
  for(unsigned int i=0;i<dominates.size();i++)
    dominates[i] = i;

  POIindex poii(poic);

  #pragma omp parallel for
  for(unsigned int i=0;i<poic.size();i++){
    if(dominates[i]!=i)                                   //Skip those already dominated
      continue;
    auto closest_n = poii.query(i);
    if(closest_n.size()==0)
      std::cerr<<"Nothing closest!"<<std::endl;
    #pragma omp critical
    for(const auto &n: closest_n){
      if(dominates[n]==n && dom_checker(poic[i],poic[n])) //Is n not already dominated? Does i dominate n?
        dominates[n]=i;                                   //Make i dominate n.second
    }
  }

  return dominates;
}

std::ofstream& PrintPOI(std::ofstream& fout, const POICollection &poic, const int i){
  fout<<i<<","
//          <<pois[i].cluster             <<","
      <<poic[i].overlaps.to_string()<<","
      <<poic[i].overlaps.count()    <<","
      <<(poic[i].pole.y*RAD_TO_DEG) <<","
      <<(poic[i].pole.x*RAD_TO_DEG) <<","
      <<(poic[i].rtheta*RAD_TO_DEG) <<","
      <<poic[i].mindist             <<","
      <<poic[i].maxdist             <<","
      <<(poic[i].avgdist/12)        <<","
      <<poic[i].edge_overlaps
      <<"\n";
  return fout;
}

std::ofstream& PrintPOICoordinates(std::ofstream& fout, const POICollection &poic, const int pn){
  IcosaXY p(poic[pn].pole, poic[pn].rtheta);
  for(unsigned int i=0;i<p.v.size();i++)
    fout<<pn                    <<","
//            <<poic[pn].cluster      <<","
        <<p.v[i].y*RAD_TO_DEG<<","
        <<p.v[i].x*RAD_TO_DEG<<","
        <<poic[pn].overlaps.test(i)
        <<"\n";
  return fout;
}

std::vector< std::vector<size_t> > FindNearbyOrientations(const POICollection &poic){
  Timer tmr;
  
  std::cerr<<"Finding nearby orientations..."<<std::endl;
  
  Timer tmr_bi;
  std::cerr<<"Building kd-tree"<<std::endl;
  POIindex poii(poic);
  std::cerr<<"Time = "<<tmr_bi.elapsed()<<std::endl;

  std::vector< std::vector<size_t> > oneighbors;
  oneighbors.reserve(poic.size());
  
  #pragma omp parallel for
  for(unsigned int i=0;i<poic.size();i++){
    const auto temp = poii.query(i);
    #pragma omp critical
    oneighbors.push_back(temp);
  }

  std::cerr<<"Finished. Time = "<<tmr_bi.elapsed()<<std::endl;

  return oneighbors;
}

TEST_CASE("POIindex"){
  POICollection poic;
  poic.emplace_back(std::bitset<12>(), Point2D(-93,45).toRadians(), 0);
  poic.emplace_back(std::bitset<12>(), Point2D(-93.1,45.1).toRadians(), 0*DEG_TO_RAD);
  poic.emplace_back(std::bitset<12>(), Point2D(-93.1,45.1).toRadians(), 72*DEG_TO_RAD);
  poic.emplace_back(std::bitset<12>(), Point2D(-93.2,45.1).toRadians(), 0*DEG_TO_RAD);
  poic.emplace_back(std::bitset<12>(), Point2D(-93.1,45.2).toRadians(), 0*DEG_TO_RAD);
  poic.emplace_back(std::bitset<12>(), Point2D(-93.2,45.2).toRadians(), 0*DEG_TO_RAD);
  poic.emplace_back(std::bitset<12>(), Point2D(-93.2,45.2).toRadians(), 36*DEG_TO_RAD);
  poic.emplace_back(std::bitset<12>(), Point2D(23,-23.2).toRadians(), 36*DEG_TO_RAD);
  auto oneighbors = FindNearbyOrientations(poic);
  CHECK(oneighbors[0].size()==5);
  std::cerr<<"oneighbors = ";
  for(const auto &x: oneighbors[7])
    std::cerr<<x<<" ";
  std::cerr<<std::endl;
  CHECK(oneighbors[7].size()==0);
}


/*
void DetermineDominants(POICollection &poic){
  Timer tmr;

  std::cerr<<"Determining dominants..."<<std::endl;

  std::cerr<<"Building POI kd-tree index..."<<std::endl;
  Timer tmr_bi;
  POIindex poii(poic);
  std::cerr<<"Finished. Time = "<<tmr_bi.elapsed()<<std::endl;

  std::cerr<<"Using tree to find nearby orientations..."<<std::endl;  


  {
    std::ofstream fout_td("test_dom");
    auto dom_checker = [](const POI &a, const POI &b){
      return a.mindist>b.mindist;
    };
    auto result = Dominants(poic, dom_checker);
    int count   = 0;
    for(unsigned int p=0;p<poic.size();p++){
      if(p!=result[p])
        continue;
      count++;
      PrintPOI(fout_td, poic, p);
    }
    std::cerr<<"Dominants found = "<<count<<std::endl;
  }

  // {
  //   auto dom_checker = [](const POI &a, const POI &b){
  //     return a.maxdist>b.maxdist;
  //   };
  //   auto result = Dominants(poic, dom_checker);
  // }

  // {
  //   auto dom_checker = [](const POI &a, const POI &b){
  //     return a.avgdist>b.avgdist;
  //   };
  //   auto result = Dominants(poic, dom_checker);
  // }

  // {
  //   auto dom_checker = [](const POI &a, const POI &b){
  //     return a.edge_overlaps<b.edge_overlaps;
  //   };
  //   auto result = Dominants(poic, dom_checker);
  // }

  std::cerr<<"Time taken = "<<tmr.elapsed()<<std::endl;
}*/

#ifdef DOCTEST_CONFIG_DISABLE

int main(){
  std::cerr<<"PRECISION = "<<PRECISION<<std::endl;

  POICollection poic;
  if(!LoadPOICollection(poic,"poic.save")){

    auto landmass = IndexedShapefile(FILE_MERC_LANDMASS,"land_polygons");

    TestWithData(landmass);

    poic = FindOrientationsOfInterest(landmass);

    DistancesToIcosaXYs(poic);

    EdgeOverlaps(landmass, poic);

    SavePOICollection(poic, "poic.save");
  }

  DetermineDominants(poic);

  std::cerr<<"Writing output..."<<std::endl;
  {
    std::ofstream fout(FILE_OUTPUT_ROT);
    fout<<"Num,Cluster,Overlaps,OverlapCount,Lat,Lon,Theta,MinDistance,MaxDistance,AvgDistance,EdgeOverlaps\n";
    for(unsigned int i=0;i<poic.size();i++)
      PrintPOI(fout,poic,i);
  }

  {
    std::ofstream fout(FILE_OUTPUT_VERT);
    fout<<"Num,Cluster,Lat,Lon,OnLand\n";
    for(unsigned int pn=0;pn<poic.size();pn++)
      PrintPOICoordinates(fout, poic, pn);
  }

  return 0;
}

#endif