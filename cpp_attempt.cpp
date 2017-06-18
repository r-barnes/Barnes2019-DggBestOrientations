//https://gis.stackexchange.com/questions/70800/how-to-read-vector-datasets-using-gdal-library
//https://math.stackexchange.com/questions/44080/approximating-geodesic-distances-on-the-sphere-by-euclidean-distances-of-a-trans

#ifndef DOCTEST_CONFIG_DISABLE
  #define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#endif

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
const double PSPACE = 10;  //km - Desired interpoint spacing

const double DEG_TO_RAD = M_PI/180.0;
const double RAD_TO_DEG = 180.0/M_PI;

//1/10th of a degree grid spacing
const double PRECISION  = 0.1;

//Neighbouring vertices of icosahedron used for generating rotations
const int NA = 0;
const int NB = 2;
const int NC = 4;

typedef std::vector< std::vector<unsigned int> > norientations_t;

template<class T>
bool LoadFromArchive(T &poic, std::string filename){
  std::ifstream os(filename, std::ios::binary);
  if(!os.good())
    return false;
  cereal::BinaryInputArchive archive( os );
  archive(poic);
  return true;
}

template<class T>
void SaveToArchive(const T &poic, std::string filename){
  std::ofstream os(filename, std::ios::binary);
  cereal::BinaryOutputArchive archive( os );
  archive(poic);
}

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

std::vector<Point2D> GenerateOrientations(const double spacing, const double radial_limit){
  std::cerr<<"Generating orientations..."<<std::endl;
  std::vector<Point2D> orientations;

  //Number of points to sample
  const int N = (int)(8*M_PI*Rearth*Rearth/std::sqrt(3)/spacing/spacing);

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

    if(temp.y<M_PI/2-radial_limit)
      continue;

    // if(!(0 <= temp.x && temp.x<=78*DEG_TO_RAD))
    //   continue;
    // if(temp.y<23*DEG_TO_RAD)
    //   continue;

    orientations.push_back(temp);
  }

  return orientations;
}

std::vector<Point2D> GenerateNearbyOrientations(const Point2D &p2d, const double angular_radius, const double spacingkm){
  auto orientations = GenerateOrientations(spacingkm,angular_radius);
  const Rotator r(Point3D(0,0,1), p2d.toXYZ(1)); //Rotates from North Pole to alternate location
  
  for(auto &x: orientations)
    x = r(x.toXYZ(1)).toLatLon();

  return orientations;
}

TEST_CASE("GenerateNearbyOrientations"){
  const auto focal            = Point2D(-93,45).toRadians();
  const double angular_radius = 2*DEG_TO_RAD;
  const auto orientations     = GenerateNearbyOrientations(focal, angular_radius, 10);

  {
    std::ofstream fout("test_nearby_orientations.csv");
    fout<<"lon,lat\n";
    for(const auto &x: orientations)
      fout<<(x.x*RAD_TO_DEG)<<","<<(x.y*RAD_TO_DEG)<<"\n";
  }

  //Check that distances to all points are within the desired angular radius of
  //the specified focal point, to within a 5% tolerance
  double maxdist = -std::numeric_limits<double>::infinity();
  for(const auto &x: orientations){
    const auto dist = GeoDistanceHaversine(focal,x);
    maxdist = std::max(maxdist,dist);
    CHECK(dist<angular_radius*6371*1.05);
  }
  std::cerr<<"Maximum nearby rotated orientation distance = "<<maxdist<<std::endl;
}

POICollection FindOrientationsOfInterest(const IndexedShapefile &landmass){
  std::cerr<<"Finding poles..."<<std::endl;
  POICollection poic;
  
  const auto orientations = GenerateOrientations(PSPACE, 90*DEG_TO_RAD);
  
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

//Maximum value returned is num_pts+1
unsigned int EdgeOverlapHelper(const IndexedShapefile &landmass, const Point2D &a, const Point2D &b, const int num_pts){
  static const GeographicLib::Geodesic geod(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f());
  unsigned int edge_overlaps = 0;
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
    edge_overlaps += PointOverlaps(temp,landmass);
  }
  return edge_overlaps;
}

void EdgeOverlaps(const IndexedShapefile &landmass, POICollection &poic){
  Timer tmr;
  std::cerr<<"Calculating edge overlaps..."<<std::endl;
  const auto   neighbors = IcosaXY().neighbors();
  const double ndist     = IcosaXY().neighborDistance()*1000; //Approximate spacing between vertices in metres
  const double spacing   = 10e3;                              //Spacing between points = 10km
  const int    num_pts   = int(std::ceil(ndist / spacing));   //The number of intervals
  std::cerr<<"Using "<<num_pts<<" with a "<<spacing<<"m spacing to cover "<<ndist<<"m inter-neighbour distance."<<std::endl;
  #pragma omp parallel for default(none) shared(poic,landmass)
  for(unsigned int pn=0;pn<poic.size();pn++){
    IcosaXY p(poic[pn].pole, poic[pn].rtheta);
    for(unsigned int n=0;n<neighbors.size();n+=2){
      const auto &a = p.v[neighbors[n]];
      const auto &b = p.v[neighbors[n+1]];
      poic[pn].edge_overlaps += EdgeOverlapHelper(landmass, a, b, num_pts);
    }
  }
  std::cout << "Time taken = " << tmr.elapsed() <<"s"<< std::endl;
}

//CHEESE
template<class T>
std::vector<unsigned int> Dominants(
  const norientations_t &orientations,
  const POICollection &poic,
  T dom_checker
){
  std::vector<size_t> dominates(poic.size());
  for(unsigned int i=0;i<dominates.size();i++)
    dominates[i] = i;

  POIindex poii(poic);

  #pragma omp parallel for default(none) schedule(static) shared(orientations,std::cerr,poic,dom_checker,dominates)
  for(unsigned int i=0;i<orientations.size();i++){
    #pragma omp critical
    for(const auto &n: orientations[i]){
      //n is already dominated
      if(dominates[n]!=n) 
        continue;
      //Don't dominate orientations with differing coverages
      if(poic[i].overlaps.count()!=poic[n].overlaps.count()) 
        continue;
      //If i doesn't dominate n, then it doesn't
      if(!dom_checker(poic[i],poic[n])) 
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

std::ofstream& PrintPOI(std::ofstream& fout, const POICollection &poic, const int i, bool header){
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
  }
  return fout;
}

std::ofstream& PrintPOICoordinates(std::ofstream& fout, const POICollection &poic, const int pn, bool header){
  if(header){
    fout<<"Num,Lat,Lon,OnLand\n";
  } else {
    IcosaXY p(poic[pn].pole, poic[pn].rtheta);
    for(unsigned int i=0;i<p.v.size();i++)
      fout<<pn                    <<","
          <<p.v[i].y*RAD_TO_DEG<<","
          <<p.v[i].x*RAD_TO_DEG<<","
          <<poic[pn].overlaps.test(i)
          <<"\n";
  }
  return fout;
}

norientations_t FindNearbyOrientations(const POICollection &poic){
  Timer tmr_bi;
  std::cerr<<"Building kd-tree"<<std::endl;
  POIindex poii(poic);
  std::cerr<<"Time = "<<tmr_bi.elapsed()<<std::endl;

  Timer tmr;
  std::cerr<<"Finding nearby orientations..."<<std::endl;

  norientations_t oneighbors(poic.size());

  #pragma omp parallel for default(none) schedule(static) shared(poic,oneighbors,poii)
  for(unsigned int i=0;i<poic.size();i++)
    oneighbors[i] = poii.query(i);

  std::cerr<<"Time taken = "<<tmr.elapsed()<<std::endl;

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
  CHECK(oneighbors[7].size()==0);
}

void PrintOrienations(
  std::string fileprefix,
  const POICollection &poic,
  const std::vector<unsigned int> &to_print
){
  std::cerr<<"Printing "<<to_print.size()<<" to '"<<fileprefix<<"'"<<std::endl;
  {
    std::ofstream fout(fileprefix+"-rot.csv");
    PrintPOI(fout, poic, 0, true);
    for(const auto &x: to_print)
      PrintPOI(fout, poic, x, false);
  }
  {
    std::ofstream fout(fileprefix+"-vert.csv");
    PrintPOICoordinates(fout, poic, 0, true);
    for(const auto &x: to_print)
      PrintPOICoordinates(fout, poic, x, false);
  }
}

void DetermineDominants(POICollection &poic, const norientations_t &norientations){
  Timer tmr;
  std::cerr<<"Determining dominants..."<<std::endl;
  {
    auto dom_checker = [](const POI &a, const POI &b){ return a.mindist>b.mindist; };
    const auto result = Dominants(norientations, poic, dom_checker);
    PrintOrienations("out_mindist", poic, result);
  }

  {
    auto dom_checker = [](const POI &a, const POI &b){ return a.maxdist>b.maxdist; };
    const auto result = Dominants(norientations, poic, dom_checker);
    PrintOrienations("out_maxdist", poic, result);
  }

  {
    auto dom_checker = [](const POI &a, const POI &b){ return a.avgdist>b.avgdist; };
    const auto result = Dominants(norientations, poic, dom_checker);
    PrintOrienations("out_avgdist", poic, result);
  }

  {
    auto dom_checker = [](const POI &a, const POI &b){ return a.edge_overlaps<b.edge_overlaps; };
    const auto result = Dominants(norientations, poic, dom_checker);
    PrintOrienations("out_edge_overlaps", poic, result);
  }

  std::cerr<<"Time taken = "<<tmr.elapsed()<<std::endl;
}



TEST_CASE("Test with data [expensive]"){
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

  SUBCASE("PointOverlaps"){
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
      CHECK(PointOverlaps(x, landmass)==true);

    //Check that all GC arcs between cities include at least some land
    for(unsigned int i=0;i<cities.size();i++)
    for(unsigned int j=0;j<cities.size();j++){
      if(i==j)
        continue;
      CHECK(EdgeOverlapHelper(landmass,cities[i],cities[j],1000)>0);
    }

    //Route from Minneapolis to Denver should be entirely on land
    {
      auto a = Point2D(-93,45).toRadians();
      auto b = Point2D(-104.9903, 39.7392).toRadians();
      CHECK(EdgeOverlapHelper(landmass,a,b,100)==101);
    }
  }
}

//Determine the number of orientations in one quadrant of the 3-space
TEST_CASE("Counting orientations [expensive]"){
  const auto orientations = GenerateOrientations(200,90*DEG_TO_RAD);

  CHECK(orientations.size()>0);

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

TEST_CASE("POIindex: Load and Save"){
  const auto a = POI(std::bitset<12>(), Point2D(-93.1,45.1).toRadians(), 1*DEG_TO_RAD);
  const auto b = POI(std::bitset<12>(), Point2D(176,-10.1).toRadians(), 7*DEG_TO_RAD);
  const auto c = POI(std::bitset<12>(), Point2D(72.4,89.3).toRadians(), 34*DEG_TO_RAD);
  const auto d = POI(std::bitset<12>(), Point2D(-103.2,-41.2).toRadians(), 98*DEG_TO_RAD);

  {
    POICollection poic;
    poic.push_back(a);
    poic.push_back(b);
    poic.push_back(c);
    poic.push_back(d);
    SaveToArchive(poic, "ztest_poic");
  }

  {
    POICollection poic;
    CHECK(LoadFromArchive(poic,"asdfasfjkwefjewifj")==false);
    CHECK(LoadFromArchive(poic,"ztest_poic")==true);
    CHECK(poic[0].pole.x==a.pole.x);
    CHECK(poic[1].pole.x==b.pole.x);
    CHECK(poic[2].pole.x==c.pole.x);
    CHECK(poic[3].pole.x==d.pole.x);
  }
}



#ifdef DOCTEST_CONFIG_DISABLE

int main(){
  std::cerr<<"PRECISION = "<<PRECISION<<std::endl;

  POICollection poic;
  if(!LoadFromArchive(poic,"poic.save")){

    auto landmass = IndexedShapefile(FILE_MERC_LANDMASS,"land_polygons");

    poic = FindOrientationsOfInterest(landmass);

    DistancesToIcosaXYs(poic);

    EdgeOverlaps(landmass, poic);

    SaveToArchive(poic, "poic.save");
  }

  norientations_t norientations;
  if(!LoadFromArchive(norientations,"norientations.save")){
    norientations = FindNearbyOrientations(poic);
    SaveToArchive(norientations, "norientations.save");
  }

  DetermineDominants(poic, norientations);

  std::cerr<<"Writing output..."<<std::endl;
  {
    std::ofstream fout(FILE_OUTPUT_ROT);
    PrintPOI(fout,poic,0,true);
    for(unsigned int i=0;i<poic.size();i++)
      PrintPOI(fout,poic,i,false);
  }

  {
    std::ofstream fout(FILE_OUTPUT_VERT);
    PrintPOICoordinates(fout, poic, 0, true);
    for(unsigned int pn=0;pn<poic.size();pn++)
      PrintPOICoordinates(fout, poic, pn, false);
  }

  return 0;
}

#endif