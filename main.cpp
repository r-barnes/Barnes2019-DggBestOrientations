//https://gis.stackexchange.com/questions/70800/how-to-read-vector-datasets-using-gdal-library
//https://math.stackexchange.com/questions/44080/approximating-geodesic-distances-on-the-sphere-by-euclidean-distances-of-a-trans

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
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicLine.hpp>
#include <GeographicLib/Constants.hpp>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <array>
#include <bitset>
#include <iomanip>
#include <cassert>
#include <fstream>
#include "doctest.h"
#include <omp.h>

#ifdef ENV_XSEDE
  const std::string FILE_WGS84_LANDMASS = "/home/rbarnes1/scratch/dgg_best/land-polygons-complete-4326/land_polygons.shp";
  const std::string FILE_OUTPUT_PREFIX  = "/home/rbarnes1/scratch/dgg_best/";
  const std::string FILE_MERC_LANDMASS  = "/home/rbarnes1/scratch/dgg_best/land-polygons-split-3857/land_polygons.shp";
#elif ENV_LAPTOP
  const std::string FILE_WGS84_LANDMASS = "data/land-polygons-complete-4326/land_polygons.shp";
  const std::string FILE_OUTPUT_PREFIX  = "/z/";
  const std::string FILE_MERC_LANDMASS  = "data/land-polygons-split-3857/land_polygons.shp";
#else
  this-is-an-error
#endif

const double Rearth = 6371; //km

const double DEG_TO_RAD = M_PI/180.0;
const double RAD_TO_DEG = 180.0/M_PI;

#ifdef FINE_RESOLUTION //Used for science
  const double COARSE_SPACING      = 10;  //km - Desired interpoint spacing for finding prospective orienations
  const double COARSE_RADIAL_LIMIT = 90*DEG_TO_RAD;
  const double COARSE_THETA_MIN    = 0;
  const double COARSE_THETA_MAX    = 72*DEG_TO_RAD;
  const double COARSE_THETA_STEP   = 0.1*DEG_TO_RAD;

  const double FINE_SPACING        = 0.1; //km - Desired interpoint spacing for zooming in on orientations of interest
  const double FINE_RADIAL_LIMIT   = COARSE_SPACING/Rearth;
  const double FINE_THETA_STEP     = 0.001*DEG_TO_RAD;
  const double FINE_THETA_INTERVAL = COARSE_THETA_STEP;
#elif COARSE_RESOLUTION //Used for profiling
  const double COARSE_SPACING      = 200;  //km - Desired interpoint spacing for finding prospective orienations
  const double COARSE_RADIAL_LIMIT = 90*DEG_TO_RAD;
  const double COARSE_THETA_MIN    = 0;
  const double COARSE_THETA_MAX    = 72*DEG_TO_RAD;
  const double COARSE_THETA_STEP   = 1.0*DEG_TO_RAD;

  const double FINE_SPACING        = 50; //km - Desired interpoint spacing for zooming in on orientations of interest
  const double FINE_RADIAL_LIMIT   = COARSE_SPACING/Rearth;
  const double FINE_THETA_STEP     = 0.5*DEG_TO_RAD;
  const double FINE_THETA_INTERVAL = COARSE_THETA_STEP;
#else
  this-is-an-error
#endif

const int    OVERLAPS_TO_BEAT    = 8; //Number of overlaps beyond (and including) which we are interested

typedef std::vector< std::vector<unsigned int> > norientations_t;

template<class T>
bool LoadFromArchive(T &poic, std::string filename){
  std::ifstream os(filename, std::ios::binary);
  if(!os.good())
    return false;
  std::cerr<<"Loading from archive '"<<filename<<"'..."<<std::endl;
  Timer tmr;
  cereal::BinaryInputArchive archive( os );
  archive(poic);
  std::cerr<<"Time = "<<tmr.elapsed()<<" s"<<std::endl;
  return true;
}

template<class T>
void SaveToArchive(const T &poic, std::string filename){
  std::cerr<<"Saving to archive '"<<filename<<"'..."<<std::endl;
  std::ofstream os(filename, std::ios::binary);
  cereal::BinaryOutputArchive archive( os );
  archive(poic);
}



PointCloud ReadPointCloudFromShapefile(std::string filename, std::string layer){
  std::cerr<<"Generating point cloud from shapefile..."<<std::endl;
  auto landmass_wgs84 = ReadShapefile(filename, layer);
  for(auto &ll: landmass_wgs84)
    ll.toRadians();

  PointCloud wgs84pc;
  std::cerr<<"Adding polygon points to PointCloud..."<<std::endl;
  for(const auto &poly: landmass_wgs84)
  for(const auto &ll: poly.exterior)
    wgs84pc.addPoint(ll.toXYZ(1));
  std::cerr<<"Building kd-tree..."<<std::endl;

  return wgs84pc;
}



//Orientations with no vertices on land are always of interest, as are overlaps
//with many orientations on land. This convenience function is used elsewhere to
//filter by overlaps
bool OverlapOfInterest(const std::bitset<OrientationWithStats::dim> &overlaps){
  return overlaps.count()==0 || overlaps.count()>=OVERLAPS_TO_BEAT;
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
auto OrientationOverlaps(const SolidXY &i2d, const IndexedShapefile &landmass){
  std::bitset<SolidXY::verts> overlaps = 0;
  for(unsigned int vi=0;vi<i2d.v.size();vi++)
    if(PointInLandmass(i2d.v[vi],landmass))
      overlaps.set(vi); 
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

  Timer tmr;
  OCollection ret;
  #pragma omp declare reduction (merge : OCollection : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
  #pragma omp parallel for default(none) schedule(static) shared(landmass) reduction(merge: ret)
  for(unsigned int o=0;o<og.size();o++){
  for(double theta=COARSE_THETA_MIN;theta<=COARSE_THETA_MAX;theta+=COARSE_THETA_STEP){
    auto ori = Orientation(og(o),theta);
    SolidXY sxy(ori);
    const auto overlaps = OrientationOverlaps(sxy, landmass);
    if(OverlapOfInterest(overlaps))
      ret.push_back(ori);
  }

  std::cout<<"Filtering: Time taken = " <<tmr.elapsed() <<"s"<<std::endl;
  std::cerr<<"Filtering: Found = "      <<ret.size()         <<std::endl;

  return ret;
}



//For a great circle connecting two points, generate sample points along the
//circle. Then count how many of these sample points fall within landmasses.
//Maximum value returned is `num_pts+1`
unsigned int GreatCircleOverlaps(const IndexedShapefile &landmass, const Point2D &a, const Point2D &b, const int num_pts){
  static const GeographicLib::Geodesic geod(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f());
  const GeographicLib::GeodesicLine line = geod.InverseLine(
    a.y*RAD_TO_DEG,
    a.x*RAD_TO_DEG,
    b.y*RAD_TO_DEG,
    b.x*RAD_TO_DEG
  );
  const double da = line.Arc() / num_pts;
  unsigned int edge_overlaps = 0;
  for(int i=0;i<=num_pts;i++) {
    Point2D temp;
    line.ArcPosition(i * da, temp.y, temp.x);
    temp.toRadians();
    edge_overlaps += PointInLandmass(temp,landmass);
  }
  return edge_overlaps;
}



//For all the great circle edges of a polyhedron, determine how many sample
//points along the circles fall within landmasses.
unsigned int OrientationEdgeOverlaps(const SolidXY &sxy, const IndexedShapefile &landmass){
  static const auto   neighbors = SolidXY().neighbors();             //Get a list of neighbouring vertices on the polyhedron
  static const double ndist     = SolidXY().neighborDistance()*1000; //Approximate spacing between vertices in metres
  static const double spacing   = 10e3;                              //Spacing between points = 10km
  static const int    num_pts   = int(std::ceil(ndist / spacing));   //The number of intervals
  unsigned int edge_overlaps = 0;
  for(unsigned int n=0;n<neighbors.size();n+=2){
    const auto &a = sxy.v[neighbors[n]];
    const auto &b = sxy.v[neighbors[n+1]];
    edge_overlaps += GreatCircleOverlaps(landmass, a, b, num_pts);
  }
  return edge_overlaps;
}



//Generate distance statistic for an orientation
OrientationWithStats OrientationStats(const Orientation &o, const PointCloud &wgs84pc, const IndexedShapefile &landmass){
  OrientationWithStats ows(o);
  SolidXY sxy(o);

  ows.overlaps = OrientationOverlaps(sxy, landmass);

  for(unsigned int i=0;i<sxy.v.size();i++){
    const auto cp = wgs84pc.queryPoint(sxy.v[i].toXYZ(1)); //Closest point
    auto llc      = cp.toLatLon();
    auto dist     = GeoDistanceFlatEarth(llc,sxy.v[i]);
    if(ows.overlaps.test(i))
      dist = -dist;
    ows.mindist = std::min(ows.mindist,dist);
    ows.maxdist = std::max(ows.maxdist,dist);
    ows.avgdist += dist;
  }
  ows.edge_overlaps = OrientationEdgeOverlaps(sxy, landmass);

  return ows;
}



OSCollection GetStatsForOrientations(const OCollection &orients, const PointCloud &wgs84pc, const IndexedShapefile &landmass){
  OSCollection osc;
  Timer tmr;
  std::cerr<<"Calculating orientation statistics..."<<std::endl;
  #pragma omp declare reduction (merge : OSCollection : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
  #pragma omp parallel for default(none) shared(wgs84pc,landmass,orients) reduction(merge:osc)
  for(unsigned int pn=0;pn<orients.size();pn++)
    osc.push_back(OrientationStats(orients[pn],wgs84pc,landmass));
  std::cout << "Time taken = " << tmr.elapsed() <<"s"<< std::endl;
  return osc;
}



//For a set of orientations and their neighbours, determine which orientations
//locally maximize a criterion specified by a function `dom_checker`.
template<class T>
std::vector<unsigned int> Dominants(
  const OSCollection &osc,
  const norientations_t &orientations,
  T dom_checker
){
  std::vector<size_t> dominates(osc.size());
  for(unsigned int i=0;i<dominates.size();i++)
    dominates[i] = i;

  #pragma omp parallel for default(none) schedule(static) shared(orientations,osc,dom_checker,dominates)
  for(unsigned int i=0;i<orientations.size();i++){
    for(const auto &n: orientations[i]){
      //n is already dominated
      if(dominates[n]!=n) 
        continue;
      //Don't dominate orientations with differing coverages
      if(osc[i].overlaps.count()!=osc[n].overlaps.count()) 
        continue;
      //If i doesn't dominate n, then it doesn't
      if(!dom_checker(osc[i],osc[n])) 
        continue;
      //Make i dominate n
      dominates[n] = i;                                 
    }
  }

  std::vector<unsigned int> ret;
  for(unsigned int i=0;i<dominates.size();i++)
    if(dominates[i]==i) //Nothing dominates *this* point
      ret.push_back(i);

  return ret;
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
        <<osc[i].overlaps.to_string()<<","
        <<osc[i].overlaps.count()    <<","
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
          <<osc[pn].overlaps.test(i)
          <<"\n";
  }
  return fout;
}



void PrintOrientations(
  std::string fileprefix,
  const OSCollection &osc,
  const std::vector<unsigned int> &do_print
){
  #pragma omp critical
  std::cerr<<"Printing "<<do_print.size()<<" to '"<<fileprefix<<"'"<<std::endl;

  {
    std::ofstream fout(FILE_OUTPUT_PREFIX + fileprefix + "-rot.csv");
    PrintPOI(fout, osc, 0, true);
    for(unsigned int dp=0;dp<do_print.size();dp++)
      PrintPOI(fout, osc, do_print[dp], false);
  }
  {
    std::ofstream fout(FILE_OUTPUT_PREFIX + fileprefix + "-vert.csv");
    PrintPOICoordinates(fout, osc, 0, true);
    for(unsigned int dp=0;dp<do_print.size();dp++)
      PrintPOICoordinates(fout, osc, do_print[dp], false);
  }
}



template<class T>
void DetermineDominantsHelper(
  const std::string      fileprefix,
  const OSCollection     &osc,
  const norientations_t  &norientations,
  T dom_checker
){
  const auto dominants = Dominants(osc, norientations, dom_checker);
  std::cerr<<"Dominants size ("<<fileprefix<<") = "<<dominants.size()<<std::endl;
  PrintOrientations(fileprefix, osc, dominants);
}



void DetermineDominants(
  OSCollection           &osc,
  const norientations_t  &norientations
){
  Timer tmr;
  std::cerr<<"Determining dominants..."<<std::endl;

  //#pragma omp parallel sections 
  {
    //#pragma omp section  
    DetermineDominantsHelper("out_min_mindist", osc, norientations,
      [](const OrientationWithStats &a, const OrientationWithStats &b){ return a.mindist<b.mindist; }
    );
    //#pragma omp section    
    DetermineDominantsHelper("out_max_mindist", osc, norientations,
      [](const OrientationWithStats &a, const OrientationWithStats &b){ return a.mindist>b.mindist; }
    );
    

    //#pragma omp section
    DetermineDominantsHelper("out_min_maxdist", osc, norientations,
      [](const OrientationWithStats &a, const OrientationWithStats &b){ return a.maxdist<b.maxdist; }
    );
    //#pragma omp section    
    DetermineDominantsHelper("out_max_maxdist", osc, norientations,
      [](const OrientationWithStats &a, const OrientationWithStats &b){ return a.maxdist>b.maxdist; }
    );

    
    //#pragma omp section
    DetermineDominantsHelper("out_min_avgdist", osc, norientations,
      [](const OrientationWithStats &a, const OrientationWithStats &b){ return a.avgdist<b.avgdist; }
    );
    //#pragma omp section    
    DetermineDominantsHelper("out_max_avgdist", osc, norientations,
      [](const OrientationWithStats &a, const OrientationWithStats &b){ return a.avgdist>b.avgdist; }
    );
    

    //#pragma omp section
    DetermineDominantsHelper("out_min_edge_overlaps", osc, norientations,
      [](const OrientationWithStats &a, const OrientationWithStats &b){ return a.edge_overlaps<b.edge_overlaps; }
    );
    //#pragma omp section    
    DetermineDominantsHelper("out_max_edge_overlaps", osc, norientations,
      [](const OrientationWithStats &a, const OrientationWithStats &b){ return a.edge_overlaps>b.edge_overlaps; }
    );
  }
  
  std::cerr<<"Time taken = "<<tmr.elapsed()<<std::endl;
}



//Determine which orientations are in the local neighbourhood of another
//orientation
template<class T>
norientations_t FindNearbyOrientations(const T &osc){
  Timer tmr_bi;
  std::cerr<<"Building kd-tree"<<std::endl;
  OrientationIndex oidx(osc);
  std::cerr<<"Time = "<<tmr_bi.elapsed()<<std::endl;

  Timer tmr;
  std::cerr<<"Finding nearby orientations..."<<std::endl;

  norientations_t oneighbors(osc.size());
  #pragma omp parallel for default(none) schedule(static) shared(osc,oneighbors,oidx)
  for(unsigned int i=0;i<osc.size();i++)
    oneighbors[i] = oidx.query(i);

  std::cerr<<"Time taken = "<<tmr.elapsed()<<std::endl;

  return oneighbors;
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

    //Check that all GC arcs between cities include at least some land
    for(unsigned int i=0;i<cities.size();i++)
    for(unsigned int j=0;j<cities.size();j++){
      if(i==j)
        continue;
      CHECK(GreatCircleOverlaps(landmass,cities[i],cities[j],1000)>0);
    }

    //Route from Minneapolis to Denver should be entirely on land
    {
      auto a = Point2D(-93,45).toRadians();
      auto b = Point2D(-104.9903, 39.7392).toRadians();
      CHECK(GreatCircleOverlaps(landmass,a,b,100)==101);
    }
  }
}

TEST_CASE("POIindex: Load and Save"){
  const auto a = Orientation(Point2D(-93.1,45.1).toRadians(),    1*DEG_TO_RAD);
  const auto b = Orientation(Point2D(176,-10.1).toRadians(),     7*DEG_TO_RAD);
  const auto c = Orientation(Point2D(72.4,89.3).toRadians(),    34*DEG_TO_RAD);
  const auto d = Orientation(Point2D(-103.2,-41.2).toRadians(), 98*DEG_TO_RAD);

  {
    OCollection oc;
    oc.push_back(a);
    oc.push_back(b);
    oc.push_back(c);
    oc.push_back(d);
    SaveToArchive(oc, "test/test_oc_save");
  }

  {
    OCollection oc;
    CHECK(LoadFromArchive(oc,"asdfasfjkwefjewifj")==false);
    CHECK(LoadFromArchive(oc,"test/test_oc_save")==true);
    CHECK(oc[0].pole.x==a.pole.x);
    CHECK(oc[1].pole.x==b.pole.x);
    CHECK(oc[2].pole.x==c.pole.x);
    CHECK(oc[3].pole.x==d.pole.x);
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
  auto oneighbors = FindNearbyOrientations(oc);
  CHECK(oneighbors[0].size()==5);
  CHECK(oneighbors[7].size()==0);
}


#ifdef DOCTEST_CONFIG_DISABLE

int main(){
  std::cout<<"Rearth              = " << Rearth              <<std::endl;
  std::cout<<"COARSE_SPACING      = " << COARSE_SPACING      <<std::endl;
  std::cout<<"COARSE_RADIAL_LIMIT = " << COARSE_RADIAL_LIMIT <<std::endl;
  std::cout<<"COARSE_THETA_MIN    = " << COARSE_THETA_MIN    <<std::endl;
  std::cout<<"COARSE_THETA_MAX    = " << COARSE_THETA_MAX    <<std::endl;
  std::cout<<"COARSE_THETA_STEP   = " << COARSE_THETA_STEP   <<std::endl;

  std::cout<<"FINE_SPACING        = " << FINE_SPACING        <<std::endl;
  std::cout<<"FINE_THETA_STEP     = " << FINE_THETA_STEP     <<std::endl;
  std::cout<<"FINE_THETA_INTERVAL = " << FINE_THETA_INTERVAL <<std::endl;
  std::cout<<"FINE_RADIAL_LIMIT   = " << FINE_RADIAL_LIMIT   <<std::endl;

  assert(!omp_get_nested());

  const auto landmass = IndexedShapefile(FILE_MERC_LANDMASS,"land_polygons");

  OCollection orients;
  if(!LoadFromArchive(orients, FILE_OUTPUT_PREFIX + "save_orients.save")){  
    orients = OrientationsFilteredByOverlaps(landmass);
    SaveToArchive(orients, FILE_OUTPUT_PREFIX + "save_orients.save");
  }

  PointCloud wgs84pc;
  if(!wgs84pc.loadFromArchive(FILE_OUTPUT_PREFIX + "save_pointcloud.save")){
    wgs84pc = ReadPointCloudFromShapefile(FILE_WGS84_LANDMASS, "land_polygons");
    wgs84pc.buildIndex();
    wgs84pc.saveToArchive(FILE_OUTPUT_PREFIX + "save_pointcloud.save");
  }

  OSCollection osc;
  if(!LoadFromArchive(osc, FILE_OUTPUT_PREFIX + "save_osc.save")){
    osc = GetStatsForOrientations(orients,wgs84pc,landmass);
    SaveToArchive(osc, FILE_OUTPUT_PREFIX + "save_osc.save");
  }

  norientations_t norientations;
  if(!LoadFromArchive(norientations, FILE_OUTPUT_PREFIX + "save_norientations.save")){
    norientations = FindNearbyOrientations(osc);
    SaveToArchive(norientations, FILE_OUTPUT_PREFIX + "save_norientations.save");
  }

  DetermineDominants(osc, norientations);

  return 0;
}

#endif