#include "POI.hpp"
#include "GeoStuff.hpp"
#include <unordered_map>
#include <iostream>
#include <fstream>
#include "doctest.h"

const double DEG_TO_RAD = M_PI/180.0;
const double RAD_TO_DEG = 180.0/M_PI;



POI::POI(const std::bitset<12> &overlaps0, const Point2D &pole0, double rtheta0){
  overlaps = overlaps0;
  pole     = pole0;
  rtheta   = rtheta0;
  _ico2d   = IcosaXY(pole,rtheta);
  _ico3d   = _ico2d.toXYZ(6371); //Radius of Earth in km
}



unsigned int POI::size() const {
  return dim;
}

TEST_CASE("POI::size"){
  POI poi;
  CHECK(poi.size()==poi.ico3d().v.size());
  CHECK(poi.size()==poi.ico2d().v.size());
}

const IcosaXYZ& POI::ico3d() const{
  return _ico3d;
}

const IcosaXY&  POI::ico2d() const{
  return _ico2d;
}




inline size_t POICollection::kdtree_get_point_count() const {
  return POI::dim*pois.size();
}

// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
inline double POICollection::kdtree_distance(const double *p1, const size_t idx_p2, size_t /*size*/) const {
  const double d0 = p1[0]-pois[idx_p2/12].ico3d().v[idx_p2%12].x;
  const double d1 = p1[1]-pois[idx_p2/12].ico3d().v[idx_p2%12].y;
  const double d2 = p1[2]-pois[idx_p2/12].ico3d().v[idx_p2%12].z;
  return d0*d0+d1*d1+d2*d2;
}

// Returns the dim'th component of the idx'th point in the class:
// Since this is inlined and the "dim" argument is typically an immediate value, the
//  "if/else's" are actually solved at compile time.
inline double POICollection::kdtree_get_pt(const size_t idx, int dim) const {
  if (dim==0)
    return pois[idx/12].ico3d().v[idx%12].x;
  else if (dim==1)
    return pois[idx/12].ico3d().v[idx%12].y;
  else
    return pois[idx/12].ico3d().v[idx%12].z;
}

// Optional bounding-box computation: return false to default to a standard bbox computation loop.
//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
template <class BBOX>
bool POICollection::kdtree_get_bbox(BBOX& /* bb */) const { return false; }

void POICollection::buildIndex() {
  if(index!=NULL)
    delete index;
  index = new my_kd_tree_t(3 /*dim*/, *this, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  index->buildIndex();
}

void POICollection::addPOI(const std::bitset<12> &overlaps, const Point2D &p, double rtheta) {
  pois.emplace_back(overlaps, p, rtheta);
}

void POICollection::addPOI(const POI &poi) {
  pois.push_back(poi);
}

POI& POICollection::operator[](unsigned int i){
  return pois[i];
}

const POI& POICollection::operator[](unsigned int i) const {
  return pois[i];
}

unsigned int POICollection::size() const {
  return pois.size();
}

std::vector<size_t> POICollection::query(const Point3D &qp) const {
  double query_pt[3] = {qp.x,qp.y,qp.z};
  const size_t num_results = 10000;
  size_t ret_index[num_results];
  double out_dist_sqr[num_results];
  nanoflann::KNNResultSet<double> resultSet(num_results);
  resultSet.init(ret_index, out_dist_sqr);
  index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));
  std::vector<size_t> temp(ret_index, ret_index+num_results);
  if(size()<num_results)
    temp.resize(size());
  return temp;
}

std::vector<size_t> POICollection::query(const unsigned int qpn) const {
  const auto &qp   = pois[qpn];
  const auto &qp2d = qp.ico2d();

  //For each 3D point of the query POI, find its nearest neighbours in 3-space
  std::vector< std::vector<size_t> > results(POI::dim);
  for(unsigned int i=0;i<POI::dim;i++)
    results[i] = query(qp.ico3d().v[i]);

  //Look at the neighbours of each point and determine how many times each
  //neighbour is seen total.
  std::unordered_map<size_t,uint8_t> counts(10000);
  for(const auto &r: results)
  for(const auto &x: r){
    if(x/12!=qpn)
      counts[x/12]++;
  }

  //Neighbours seen near all 12 points represent an actual neighbouring
  //configuration, rather than just one that is nearby. Let us remove those
  //configurations which do not fit this criterion.
  for(auto it = counts.begin(); it!=counts.end(); )
    if(it->second!=12)
      it = counts.erase(it);
    else
      ++it;

  //For those that remain, sum their distances to the query point
  std::unordered_map<size_t,double> distances(10000);
  for(unsigned int i=0;i<results.size();i++)
  for(const auto &x: results[i]){
    if(counts.count(x/12)==0)
      continue;
    auto dist = GeoDistanceHaversine(pois[x/12].ico2d().v[x%12], qp2d.v[i]);
    if(dist>100)
      distances[x/12] += std::numeric_limits<double>::infinity();
    else
      distances[x/12] += dist;
  }

  //Put all neighbours which are close enough into the vector
  std::vector<size_t> closest_n;
  for(const auto &nd: distances){
    if(nd.second/12<=200)
      closest_n.emplace_back(nd.first);
  }
  std::sort(closest_n.begin(),closest_n.end(), [&](const size_t a, const size_t b){return distances[a]<distances[b];});

  if(closest_n.size()>10)
    closest_n.resize(10);

  return closest_n;
}

void POICollection::save(std::string filename) const {
  std::ofstream os(filename, std::ios::binary);
  cereal::BinaryOutputArchive archive( os );
  archive(*this);
}

bool POICollection::load(std::string filename) {
  std::ifstream os(filename, std::ios::binary);
  if(!os.good())
    return false;
  cereal::BinaryInputArchive archive( os );
  archive(*this);
  return true;
}

TEST_CASE("POICollection: Load and Save"){
  const auto a = POI(std::bitset<12>(), Point2D(-93.1,45.1).toRadians(), 1*DEG_TO_RAD);
  const auto b = POI(std::bitset<12>(), Point2D(176,-10.1).toRadians(), 7*DEG_TO_RAD);
  const auto c = POI(std::bitset<12>(), Point2D(72.4,89.3).toRadians(), 34*DEG_TO_RAD);
  const auto d = POI(std::bitset<12>(), Point2D(-103.2,-41.2).toRadians(), 98*DEG_TO_RAD);

  {
    POICollection poic;
    poic.addPOI(a);
    poic.addPOI(b);
    poic.addPOI(c);
    poic.addPOI(d);
    poic.save("ztest_poic");
  }

  {
    POICollection poic;
    CHECK(poic.load("asdfasfjkwefjewifj")==false);
    CHECK(poic.load("ztest_poic")==true);
    CHECK(poic[0].pole.x==a.pole.x);
    CHECK(poic[1].pole.x==b.pole.x);
    CHECK(poic[2].pole.x==c.pole.x);
    CHECK(poic[3].pole.x==d.pole.x);
    CHECK(poic[0].ico3d().v[3].x==a.ico3d().v[3].x);
    CHECK(poic[0].ico2d().v[5].y==a.ico2d().v[5].y);
  }
}



TEST_CASE("POICollection"){
  POICollection poic;
  poic.addPOI(std::bitset<12>(), Point2D(-93,45).toRadians(), 0);

  SUBCASE("Indexing"){
    CHECK(&poic[0] == &poic.pois[0]);
  }

  SUBCASE("No result"){
    poic.buildIndex();
    CHECK(poic.size()==1);
    auto result = poic.query(0);
    CHECK(result.size()==0);
  }

  SUBCASE("One result"){
    poic.addPOI(std::bitset<12>(), Point2D(-93,45).toRadians(), 72.0*DEG_TO_RAD);
    poic.buildIndex();
    CHECK(poic.size()==2);
    auto result = poic.query(0);
    CHECK(result[0]==1);
  }

  SUBCASE("One result: build index twice"){
    poic.addPOI(std::bitset<12>(), Point2D(-93,45).toRadians(), 72.0*DEG_TO_RAD);
    poic.buildIndex();
    poic.buildIndex();
    CHECK(poic.size()==2);
    auto result = poic.query(0);
    CHECK(result[0]==1);
  }

  SUBCASE("Check more"){
    poic.addPOI(std::bitset<12>(), Point2D(-93.1,45.1).toRadians(), 0*DEG_TO_RAD);
    poic.addPOI(std::bitset<12>(), Point2D(-93.2,45.1).toRadians(), 0*DEG_TO_RAD);
    poic.addPOI(std::bitset<12>(), Point2D(-93.1,45.2).toRadians(), 0*DEG_TO_RAD);
    poic.addPOI(std::bitset<12>(), Point2D(-93.2,45.2).toRadians(), 0*DEG_TO_RAD);
    poic.addPOI(std::bitset<12>(), Point2D(-93.2,45.2).toRadians(), 36*DEG_TO_RAD);
    poic.addPOI(std::bitset<12>(), Point2D(23,-23.2).toRadians(), 36*DEG_TO_RAD);
    poic.buildIndex();
    CHECK(poic.size()==7);
    auto result = poic.query(0);
    CHECK(result.size()==4);
    CHECK(result[0]==1);
  }
}
