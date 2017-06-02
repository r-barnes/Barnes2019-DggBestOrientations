//https://gis.stackexchange.com/questions/70800/how-to-read-vector-datasets-using-gdal-library
//https://math.stackexchange.com/questions/44080/approximating-geodesic-distances-on-the-sphere-by-euclidean-distances-of-a-trans
#include "Polygon.hpp"
#include "SpIndex.hpp"
#include "PointCloud.hpp"
#include "Icosa.hpp"
#include "GeoStuff.hpp"
#include "Point.hpp"
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
#include <bitset>

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

void AddPolygonToSpIndex(const Polygon &poly, SpIndex &sp, const int id){
  const int xmin = poly.minX();
  const int ymin = poly.minY();
  const int xmax = poly.maxX();
  const int ymax = poly.maxY();
  sp.addBoxDeferred(xmin,ymin,xmax,ymax,id);
}

void TestWithData(const Polygons &landmass_merc, const SpIndex &sp){
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
}

void Test(){
  std::cerr<<"Running tests..."<<std::endl;

  {
    std::cerr<<"Ico 2D coords in XYZ"<<std::endl;
    IcosaXY p;
    p.toXYZ().print();
  }

  {
    std::cerr<<"IcosaXY"<<std::endl;
    IcosaXY p;
    p.rotate(22.8*DEG_TO_RAD,3.6*DEG_TO_RAD,45.6*DEG_TO_RAD);
    p.print();
    p.toMercator();
    p.print();
  }

  {
    std::cerr<<"IcosoXY Neighbors:"<<std::endl;
    IcosaXY p;
    const auto n = p.neighbors();
    for(unsigned int i=0;i<n.size();i+=2)
      std::cerr<<n[i]<<"-"<<n[i+1]<<std::endl;
    //Ensure neighbours are actually all neighbours
    const auto dist = GeoDistanceHaversine(p.v[NA],p.v[NB]);
    assert(std::abs(GeoDistanceHaversine(p.v[NA],p.v[NC])-dist)<1e-6);
    assert(std::abs(GeoDistanceHaversine(p.v[NB],p.v[NC])-dist)<1e-6);
  }

  SpIndex sp;
  int id=0;
  for(double y=0;y<1000;y+=100)
  for(double x=0;x<1000;x+=100)
    sp.addBox(x,y,x+100,y+100,id++);

  assert(sp.queryPoint(Point2D(350,350))==33);
  assert(sp.queryPoint(Point2D(750,550))==57);

  Polygon p;
  p.exterior.emplace_back(1200,1200);
  p.exterior.emplace_back(1200,1300);
  p.exterior.emplace_back(1300,1300);
  p.exterior.emplace_back(1300,1200);

  AddPolygonToSpIndex(p, sp, 347);
  sp.buildIndex();

  std::cerr<<sp.queryPoint(Point2D(1250,1270))<<std::endl;

  assert(sp.queryPoint(Point2D(1250,1270))==347);

  assert(p.containsPoint(Point2D(1250,1270)));
  assert(p.containsPoint(Point2D(1243,1222)));
  assert(!p.containsPoint(Point2D(1194,1222)));

  {
    Point2D ll(-93*DEG_TO_RAD,45*DEG_TO_RAD);
    ll = WGS84toEPSG3857(ll);
    std::cout<<"(-93,45) = ("<<std::fixed<<ll.x<<","<<std::fixed<<ll.y<<")"<<std::endl;
    assert(std::abs(ll.x-(-10352712.6438))<1e-4);
    assert(std::abs(ll.y-5621521.48619)<1e-4);
  }

  //TODO
  {
    //Describes a regular icosahedron with edges of length 2
    IcosaXYZ ico;
    ico.rotateTo(Point3D(0,0,1));
    auto ixy = ico.toLatLon();
    //std::cerr<<"Before rotation"<<std::endl;
    //ixy.print();
    //ixy.rotate(-26.565051*DEG_TO_RAD,-90*DEG_TO_RAD,0);
    std::cerr<<"After rotation"<<std::endl;
    ixy.print();
  }

  {
    const IcosaXYZ ico;
    const Point3D A = ico.v[NA];
    const Point3D AB(
      ico.v[NB].x-ico.v[NA].x,
      ico.v[NB].y-ico.v[NA].y,
      ico.v[NB].z-ico.v[NA].z
    );
    const Point3D AC(
      ico.v[NC].x-ico.v[NA].x,
      ico.v[NC].y-ico.v[NA].y,
      ico.v[NC].z-ico.v[NA].z
    );
    std::ofstream fout("/z/dgbp_coverage_test");
    for(int rn=0;rn<=50;rn++)
    for(int sn=0;sn<=50;sn++){
      double r = rn/(double)50;
      double s = sn/(double)50;
      if(r+s>=1.1)
        continue;
      Point3D orient_to(
        A.x + r*AB.x + s*AC.x,
        A.y + r*AB.y + s*AC.y,
        A.z + r*AB.z + s*AC.z
      );
      for(int16_t rtheta=0; rtheta<72; rtheta+=1){
        IcosaXY p;
        p = p.rotateTheta(rtheta*DEG_TO_RAD).toXYZ().rotateTo(orient_to).toLatLon();
        //fout<<(p.v[0].y*RAD_TO_DEG)<<" "<<(p.v[0].x*RAD_TO_DEG)<<"\n";
        for(const auto &i: p.v)
          fout<<(i.y*RAD_TO_DEG)<<" "<<(i.x*RAD_TO_DEG)<<"\n";
      }
    }
  }

  std::cerr<<"Passed"<<std::endl;
}

class POI {
 public:
  std::bitset<12> overlaps      = 0;
  int16_t         rlat          = 0;
  int16_t         rlon          = 0;
  int16_t         rtheta        = 0;
  double          mindist       = std::numeric_limits<double>::infinity();
  double          maxdist       = -std::numeric_limits<double>::infinity();
  double          avgdist       = 0;
  uint16_t        cluster       = 0;
  int             edge_overlaps = 0;
  POI(std::bitset<12> overlaps,double rlat,double rlon,double rtheta){
    this->overlaps = overlaps;
    this->rlat     = rlat;
    this->rlon     = rlon;
    this->rtheta   = rtheta;
  }
};

std::vector<struct POI> FindOrientationsOfInterest(
  const Polygons &landmass_merc,
  const SpIndex &sp
){
  std::cerr<<"Finding poles..."<<std::endl;
  std::vector<struct POI> pois;
  long count=0;

  const IcosaXYZ ico;
  const Point3D A = ico.v[NA];
  const Point3D AB(
    ico.v[NB].x-ico.v[NA].x,
    ico.v[NB].y-ico.v[NA].y,
    ico.v[NB].z-ico.v[NA].z
  );
  const Point3D AC(
    ico.v[NC].x-ico.v[NA].x,
    ico.v[NC].y-ico.v[NA].y,
    ico.v[NC].z-ico.v[NA].z
  );

  Timer tmr;
  //#pragma omp parallel for default(none) schedule(static) shared(pois,std::cerr,sp,landmass_merc) reduction(+:count)
  for(int rn=0;rn<=100;rn++)
  for(int sn=0;sn<=100;sn++){
    int r = rn;
    int s = sn;
    if(r+s>=1){
      r = 1-r;
      s = 1-s;
    }
    Point3D orient_to(
      A.x + r/100.0 * AB.x + s/100.0 * AC.x,
      A.y + r/100.0 * AB.y + s/100.0 * AC.y,
      A.z + r/100.0 * AB.z + s/100.0 * AC.z
    );
    for(int16_t rtheta=0; rtheta<(int)(72.0*DIV); rtheta+=(int)(PRECISION*DIV)){
      count++;
      IcosaXY p;
      p = p.rotateTheta(rtheta).toXYZ().rotateTo(orient_to).toLatLon();
      std::bitset<12> overlaps = 0;
      for(unsigned int i=0;i<p.v.size();i++)
        if(PointOverlaps(p.v[i],landmass_merc,sp))
          overlaps.set(i);
      if(overlaps==0 || overlaps.count()>=8){
        #pragma omp critical
        pois.emplace_back(overlaps,r,s,rtheta);
      }
    }
  }

  double t = tmr.elapsed();
  std::cout << "Time taken = " << t <<"s"<< std::endl;

  std::cerr<<"Checked = "<<count<<std::endl;
  std::cerr<<"Found "<<pois.size()<<" poles of interest."<<std::endl;

  return pois;
}

void DistancesToIcosaXYs(std::vector<struct POI> &pois){
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
  #pragma omp parallel for default(none) shared(pois,pc)
  for(unsigned int pn=0;pn<pois.size();pn++){
    IcosaXY p(pois[pn].rlat/DIV*DEG_TO_RAD, pois[pn].rlon/DIV*DEG_TO_RAD, pois[pn].rtheta/DIV*DEG_TO_RAD);

    for(unsigned int i=0;i<p.v.size();i++){
      const auto cp = pc.queryPoint(p.v[i].toXYZ(1)); //Closest point
      auto llc      = cp.toLatLon();
      auto dist     = GeoDistanceFlatEarth(llc,p.v[i]);
      if(pois[pn].overlaps.test(i))
        dist = -dist;
      pois[pn].mindist = std::min(pois[pn].mindist,dist);
      pois[pn].maxdist = std::max(pois[pn].maxdist,dist);
      pois[pn].avgdist += dist;
    }
  }

  std::cout << "Time taken = " << tmr.elapsed() <<"s"<< std::endl;
}

void EdgeOverlaps(
  const Polygons &landmass_merc,
  const SpIndex &sp,
  std::vector<POI> &pois
){
  Timer tmr;
  std::cerr<<"Calculating edge overlaps..."<<std::endl;
  const GeographicLib::Geodesic geod(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f());
  const auto neighbors = IcosaXY().neighbors();
  const double ndist   = IcosaXY().neighborDistance()*1000;  //Approximate spacing between vertices in metres
  const double spacing = 10e3;                            //Spacing between points = 10km
  const int    num_pts = int(std::ceil(ndist / spacing)); //The number of intervals
  #pragma omp parallel for default(none) shared(pois,landmass_merc,sp)
  for(unsigned int pn=0;pn<pois.size();pn++){
    IcosaXY p(pois[pn].rlat/DIV*DEG_TO_RAD, pois[pn].rlon/DIV*DEG_TO_RAD, pois[pn].rtheta/DIV*DEG_TO_RAD);
    for(unsigned int n=0;n<neighbors.size();n++){
      const GeographicLib::GeodesicLine line = geod.InverseLine(
        p.v[n].y*RAD_TO_DEG,
        p.v[n].x*RAD_TO_DEG,
        p.v[n+1].y*RAD_TO_DEG,
        p.v[n+1].x*RAD_TO_DEG
      );
      const double da = line.Arc() / num_pts;
      for(int i=0;i<=num_pts;i++) {
        Point2D temp;
        line.ArcPosition(i * da, temp.y, temp.x);
        temp.toRadians();
        pois[pn].edge_overlaps += PointOverlaps(temp,landmass_merc,sp);
      }
    }
  }
  std::cout << "Time taken = " << tmr.elapsed() <<"s"<< std::endl;
}

int main(int argc, char **argv){
  Test();

  std::cerr<<"PRECISION = "<<PRECISION<<std::endl;

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

  auto pois = FindOrientationsOfInterest(landmass_merc,sp);

  DistancesToIcosaXYs(pois);

  EdgeOverlaps(landmass_merc, sp, pois);

  std::cerr<<"Sorting poles of interest..."<<std::endl;
  std::sort(pois.begin(),pois.end(), [](const POI &a, const POI &b){
    return a.mindist>b.mindist;
  });

  std::cerr<<"Writing output..."<<std::endl;
  {
    std::ofstream fout(FILE_OUTPUT_ROT);
    fout<<"Num,Cluster,Overlaps,OverlapCount,Lat,Lon,Theta,MinDistance,MaxDistance,AvgDistance,EdgeOverlaps\n";
    for(unsigned int i=0;i<pois.size();i++)
      fout<<i<<","
          <<pois[i].cluster             <<","
          <<pois[i].overlaps.to_string()<<","
          <<pois[i].overlaps.count()    <<","
          <<pois[i].rlat                <<","
          <<pois[i].rlon                <<","
          <<pois[i].rtheta              <<","
          <<pois[i].mindist             <<","
          <<pois[i].maxdist             <<","
          <<(pois[i].avgdist/12)        <<","
          <<pois[i].edge_overlaps
          <<"\n";
  }

  {
    std::ofstream fout(FILE_OUTPUT_VERT);
    fout<<"Num,Cluster,Lat,Lon,OnLand\n";
    for(unsigned int pn=0;pn<pois.size();pn++){
      IcosaXY p(pois[pn].rlat/DIV*DEG_TO_RAD,pois[pn].rlon/DIV*DEG_TO_RAD,pois[pn].rtheta/DIV*DEG_TO_RAD);
      for(unsigned int i=0;i<p.v.size();i++)
        fout<<pn                    <<","
            <<pois[pn].cluster      <<","
            <<p.v[i].y*RAD_TO_DEG<<","
            <<p.v[i].x*RAD_TO_DEG<<","
            <<pois[pn].overlaps.test(i)
            <<"\n";
    }
  }

  return 0;
}
