//https://gis.stackexchange.com/questions/70800/how-to-read-vector-datasets-using-gdal-library
//https://math.stackexchange.com/questions/44080/approximating-geodesic-distances-on-the-sphere-by-euclidean-distances-of-a-trans
#include "Polygon.hpp"
#include "SpIndex.hpp"
#include "PointCloud.hpp"
#include "Icosa.hpp"
#include "GeoStuff.hpp"
#include "Point.hpp"
#include "POI.hpp"
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicLine.hpp>
#include <GeographicLib/Constants.hpp>
#include <ogrsf_frmts.h>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <array>
#include <iomanip>
#include <cassert>
#include <fstream>
#include <chrono>

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

class Timer {
 private:
  typedef std::chrono::high_resolution_clock clock_;
  typedef std::chrono::duration<double, std::ratio<1> > second_;
  std::chrono::time_point<clock_> beg_;
 public:
  Timer() : beg_(clock_::now()) {}
  void reset() { beg_ = clock_::now(); }
  double elapsed() const { 
    return std::chrono::duration_cast<second_> (clock_::now() - beg_).count(); 
  }
};

void ReadShapefile(std::string filename, std::string layername, Polygons &geometries){
  GDALAllRegister();
  GDALDataset *poDS;
  poDS = (GDALDataset*) GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL);
  if( poDS == NULL ){
    std::cerr<<"Failed to open '"<<filename<<"'!"<<std::endl;
    exit( 1 );
  }

  OGRLayer *poLayer;
  poLayer = poDS->GetLayerByName(layername.c_str());
  OGRFeature *poFeature;
  poLayer->ResetReading();
  while( (poFeature = poLayer->GetNextFeature()) != NULL ){
    /*OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
    int iField;
    for( iField = 0; iField < poFDefn->GetFieldCount(); iField++ ){
      OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn( iField );
      if( poFieldDefn->GetType() == OFTInteger )
          printf( "%d,", poFeature->GetFieldAsInteger( iField ) );
      else if( poFieldDefn->GetType() == OFTInteger64 )
          printf( CPL_FRMT_GIB ",", poFeature->GetFieldAsInteger64( iField ) );
      else if( poFieldDefn->GetType() == OFTReal )
          printf( "%.3f,", poFeature->GetFieldAsDouble(iField) );
      else if( poFieldDefn->GetType() == OFTString )
          printf( "%s,", poFeature->GetFieldAsString(iField) );
      else
          printf( "%s,", poFeature->GetFieldAsString(iField) );
    }*/
    OGRGeometry *poGeometry;
    poGeometry = poFeature->GetGeometryRef();
    if (poGeometry==NULL){
      //Pass
    } else if( wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon ){
      OGRPolygon *poly = (OGRPolygon *) poGeometry;
      auto extring = poly->getExteriorRing();
      //Ignore interior rings for now: they're probably lakes
      geometries.emplace_back();
      for(int i=0;i<extring->getNumPoints();i++)
        geometries.back().exterior.emplace_back(extring->getX(i),extring->getY(i));
    } else {
      std::cerr<<"Unrecognised geometry of type: "<<wkbFlatten(poGeometry->getGeometryType())<<std::endl;
    }
    OGRFeature::DestroyFeature( poFeature );
  }
  GDALClose( poDS );
}

bool PointOverlaps(
  const Point2D &ll,
  const Polygons &landmass_merc,
  const SpIndex &sp
){
  if(ll.y>83.7*DEG_TO_RAD) //The island "83-42" as at 83.7N so anything north of this is on water
    return false;
  if(ll.y<-80*DEG_TO_RAD)  //The southmost extent of water is ~79.5S, so anything south of this is on land
    return true;
  auto xy = WGS84toEPSG3857(ll);
  const auto pid = sp.queryPoint(xy);
  if(pid==-1)
    return false;
  if(landmass_merc.at(pid).containsPoint(xy))
    return true;
  return false;
}

void TestWithData(const Polygons &landmass_merc, const SpIndex &sp){
  std::cerr<<"Running data-based tests..."<<std::endl;
  {
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
      if(PointOverlaps(v,landmass_merc,sp))
        ocount++;

    std::cerr<<"Fuller count: "<<ocount<<std::endl;
    assert(ocount==0);
  }

  std::cerr<<"Passed"<<std::endl;
}

void Test(){
  std::cerr<<"Running tests..."<<std::endl;

  {
    const GeographicLib::Geodesic geod(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f());
    IcosaXY p;
    const auto   neighbors = p.neighbors();
    const double ndist     = IcosaXY().neighborDistance()*1000;  //Approximate spacing between vertices in metres
    const double spacing   = 10e3;                               //Spacing between points = 10km
    const int    num_pts   = int(std::ceil(ndist / spacing));    //The number of intervals

    std::cerr<<"GC arcs: "<<(neighbors.size()/2)<<std::endl;

    std::ofstream fout_gc("test_gc_points");
    const auto dist = GeoDistanceHaversine(p.v[NA],p.v[NB]);
    //std::cerr<<"GC dist nominal = "<<dist<<std::endl;
    for(unsigned int n=0;n<neighbors.size();n+=2){
      const auto &a = p.v[neighbors[n]];
      const auto &b = p.v[neighbors[n+1]];
      //std::cerr<<"GC dist diff = "<<std::abs(GeoDistanceHaversine(a,b)-dist)<<std::endl;
      assert(std::abs(GeoDistanceHaversine(a,b)-dist)<1);
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
        fout_gc<<temp.y<<" "<<temp.x<<"\n";
      }
    }
  }

  {
    POICollection poic;
    poic.addPOI(std::bitset<12>(), Point2D(-93,45).toRadians(), 0);
    poic.buildIndex();
    assert(poic.size()==1);
    auto result = poic.query(0);
    assert(result.size()==0);
    std::cerr<<"x-val: "<<poic[0].ico3d.v[0].x<<std::endl;
    std::cerr<<"y-val: "<<poic[0].ico3d.v[0].y<<std::endl;
    std::cerr<<"z-val: "<<poic[0].ico3d.v[0].z<<std::endl;
  }

  {
    POICollection poic;
    poic.addPOI(std::bitset<12>(), Point2D(-93,45).toRadians(), 0);
    poic.addPOI(std::bitset<12>(), Point2D(-93,45).toRadians(), 72.0*DEG_TO_RAD);
    poic.buildIndex();
    assert(poic.size()==2);
    auto result = poic.query(0);
    assert(result[0]==1);
  }

  std::cerr<<"Passed"<<std::endl;
}

POICollection FindOrientationsOfInterest(
  const Polygons &landmass_merc,
  const SpIndex &sp
){
  std::cerr<<"Finding poles..."<<std::endl;
  POICollection poic;
  std::vector<Point2D> orientations;
  long count=0;

  //Number of points to sample
  const int N = (int)(8*M_PI*Rearth*Rearth/std::sqrt(3)/pspace/pspace);

  //Generate orientations
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


  Timer tmr;
  #pragma omp parallel for default(none) schedule(static) shared(orientations,poic,std::cerr,sp,landmass_merc) reduction(+:count)
  for(unsigned int i=0;i<orientations.size();i++)
  for(double rtheta=0;rtheta<72.01*DEG_TO_RAD;rtheta+=PRECISION*DEG_TO_RAD){
    count++;
    IcosaXY p(orientations[i],rtheta);

    std::bitset<12> overlaps = 0;
    for(unsigned int i=0;i<p.v.size();i++)
      if(PointOverlaps(p.v[i],landmass_merc,sp))
        overlaps.set(i);
    if(overlaps==0 || overlaps.count()>=8){
      #pragma omp critical
      poic.addPOI(overlaps,orientations[i],rtheta);
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
  Polygons landmass_wgs84;
  ReadShapefile(FILE_WGS84_LANDMASS, "land_polygons", landmass_wgs84);

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

void EdgeOverlaps(
  const Polygons &landmass_merc,
  const SpIndex &sp,
  POICollection &poic
){
  Timer tmr;
  std::cerr<<"Calculating edge overlaps..."<<std::endl;
  const GeographicLib::Geodesic geod(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f());
  const auto neighbors = IcosaXY().neighbors();
  const double ndist   = IcosaXY().neighborDistance()*1000;  //Approximate spacing between vertices in metres
  const double spacing = 10e3;                            //Spacing between points = 10km
  const int    num_pts = int(std::ceil(ndist / spacing)); //The number of intervals
  #pragma omp parallel for default(none) shared(poic,landmass_merc,sp)
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
        poic[pn].edge_overlaps += PointOverlaps(temp,landmass_merc,sp);
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

  for(unsigned int i=0;i<poic.size();i++){
    if(dominates[i]!=i)                                   //Skip those already dominated
      continue;
    auto closest_n = poic.query(i);
    if(closest_n.size()==0)
      std::cerr<<"Nothing closest!"<<std::endl;
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

void DetermineDominants(POICollection &poic){
  Timer tmr;

  std::cerr<<"Determining dominants..."<<std::endl;
  std::cerr<<"Building POI kd-tree index..."<<std::endl;

  {
    Timer tmr;
    poic.buildIndex();
    std::cerr<<"Finished. Time = "<<tmr.elapsed()<<std::endl;
  }

  std::cerr<<"Using tree to search for dominants..."<<std::endl;
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
}

int main(int argc, char **argv){
  Test();

  std::cerr<<"PRECISION = "<<PRECISION<<std::endl;

  POICollection poic;
  if(!poic.load("poic.save")){

    Polygons landmass_merc;
    {
      std::cerr<<"Reading Mercator split shapefile..."<<std::endl;
      Timer tmr;
      ReadShapefile(FILE_MERC_LANDMASS, "land_polygons", landmass_merc);
      std::cerr<<"Read "<<landmass_merc.size()<<" polygons."<<std::endl;
      std::cerr<<"Time taken = "<<tmr.elapsed()<<std::endl;
    }

    SpIndex sp;
    {
      std::cerr<<"Building index..."<<std::endl;
      Timer tmr;
      for(int64_t i=0;(unsigned)i<landmass_merc.size();i++)
        AddPolygonToSpIndex(landmass_merc[i], sp, i);
      sp.buildIndex();
      std::cerr<<"Index built."<<std::endl;
      std::cerr<<"Time taken = "<<tmr.elapsed()<<std::endl;
    }

    TestWithData(landmass_merc,sp);

    poic = FindOrientationsOfInterest(landmass_merc,sp);

    DistancesToIcosaXYs(poic);

    EdgeOverlaps(landmass_merc, sp, poic);

    poic.save("poic.save");
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
