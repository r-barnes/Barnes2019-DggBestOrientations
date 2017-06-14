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
}

POIindex::POIindex(const std::vector<POI> &pois){
  //Cannot be parallelized, otherwise the order of the points might get mixed up
  for(unsigned int pi=0;pi<pois.size();pi++){
    const IcosaXY  p2d = IcosaXY(pois[pi].pole,pois[pi].rtheta);
    const IcosaXYZ p3d = p2d.toXYZ(6371);  //Radius of the Earth
    for(unsigned int vi=0;vi<p3d.v.size();vi++){
      //Choose one quadrant of 3-space. Will have 2-3 members
      if(p3d.v[vi].z>=0 && p3d.v[vi].y>=0){ 
        p3ds.push_back(p3d.v[vi]);
        p2ds.push_back(p2d.v[vi]);
        pidx.push_back(pi);
      }
    }
  }

  index.reset(new my_kd_tree_t(3 /*dim*/, *this, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) ));
  index->buildIndex();
}

POIindex::~POIindex(){
  index->freeIndex();
}

inline size_t POIindex::kdtree_get_point_count() const {
  return p3ds.size();
}

// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
inline double POIindex::kdtree_distance(const double *p1, const size_t idx_p2, size_t /*size*/) const {
  const double d0 = p1[0]-p3ds[idx_p2].x;
  const double d1 = p1[1]-p3ds[idx_p2].y;
  const double d2 = p1[2]-p3ds[idx_p2].z;
  return d0*d0+d1*d1+d2*d2;
}

// Returns the dim'th component of the idx'th point in the class:
// Since this is inlined and the "dim" argument is typically an immediate value, the
//  "if/else's" are actually solved at compile time.
inline double POIindex::kdtree_get_pt(const size_t idx, int dim) const {
  if (dim==0)
    return p3ds[idx].x;
  else if (dim==1)
    return p3ds[idx].y;
  else
    return p3ds[idx].z;
}

// Optional bounding-box computation: return false to default to a standard bbox computation loop.
//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
template <class BBOX>
bool POIindex::kdtree_get_bbox(BBOX& /* bb */) const { return false; }

std::vector<unsigned int> POIindex::query(const Point3D &qp) const {
  double query_pt[3] = {qp.x,qp.y,qp.z};

  //const size_t num_results = 10000;
  //size_t ret_index[num_results];
  //double out_dist_sqr[num_results];
  //nanoflann::KNNResultSet<double> resultSet(num_results);
  //resultSet.init(ret_index, out_dist_sqr);
  //index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));
  //std::vector<size_t> temp(ret_index, ret_index+num_results);
  //if(size()<num_results)
  //  temp.resize(size());

  std::vector<std::pair<size_t, double> > matches;
  nanoflann::SearchParams params;

  index->radiusSearch(query_pt, 100*100, matches, params);
  std::sort(matches.begin(),matches.end(),[&](const std::pair<size_t, double> &a, const std::pair<size_t, double> &b){return a.second<b.second;});
  //std::cerr<<"Found "<<matches.size()<<" in radius."<<std::endl;
  //std::cerr<<"Smallest = "<<matches.front().second<<", Largest="<<matches.back().second<<std::endl;

  std::vector<unsigned int> temp(matches.size());
  for(const auto &m: matches)
    temp.emplace_back(m.first);

  return temp;
}

std::vector<unsigned int> POIindex::query(const unsigned int qpn) const {
  //Locate the first point corresponding to qpn, we'll iterate forward to
  //identify the rest
  const size_t fs = std::lower_bound(pidx.begin(), pidx.end(), qpn)-pidx.begin();

  //For each 3D point of the query POI, find its nearest neighbours in 3-space
  std::vector< std::vector<unsigned int> > results;
  results.reserve(3); //Can have 2-3 neighbours
  //Iterate through all of the points associated with qpn that we have in the
  //index, finding their nearest neighbours
  for(unsigned int fi=fs;pidx[fi]==qpn;fi++)
    results.push_back(query(p3ds[fi]));

  //Look at the neighbours of each point and determine how many times each
  //neighbour is seen total. Ignore qpn itself.
  std::unordered_map<size_t,uint8_t> counts(10000);
  for(const auto &r: results)
  for(const auto &x: r){
    if(pidx[x]!=qpn)
      counts[pidx[x]]++;
  }

  //Neighbours seen near two points represent an actual neighbouring
  //configuration, rather than just one that is nearby. Let us remove those
  //configurations which do not fit this criterion.
  for(auto it = counts.begin(); it!=counts.end(); )
    if(it->second<2)
      it = counts.erase(it);
    else
      ++it;

  //For those that remain, sum their distances to the query point
  std::unordered_map<size_t,double> distances(10000);
  //For each vertex of qpn that falls in the search zone, get the distances from
  //that vertex to all of its nearest neighbours. However, ignore those
  //neighbours which did not also appear near at least one other vertex of qpn.
  for(unsigned int i=0;i<results.size();i++)
  for(const auto &x: results[i]){
    if(counts.count(pidx[x])==0)
      continue;
    const auto dist = GeoDistanceHaversine(p2ds[x], p2ds[fs+i]);
    if(dist>100) //Is neighbor too distant? If so, trash the whole orientation
      distances[pidx[x]] += std::numeric_limits<double>::infinity();
    else
      distances[pidx[x]] += dist;
  }

  //Put all neighbours which are close enough into the vector
  std::vector<unsigned int> closest_n;
  for(const auto &nd: distances){
    if(nd.second<std::numeric_limits<double>::infinity())
      closest_n.emplace_back(nd.first);
  }
  std::sort(closest_n.begin(), closest_n.end(), [&](const size_t a, const size_t b){return distances[a]<distances[b];});

  // std::cerr<<"Closest (qpn="<<qpn<<"): "<<std::endl;
  // for(const auto &nd: closest_n){
  //   std::cerr<<"\t"<<nd<<" "<<distances[nd]<<std::endl;
  // }

  if(closest_n.size()>10)
    closest_n.resize(10);

  return closest_n;
}



bool LoadPOICollection(POICollection &poic, std::string filename){
  std::ifstream os(filename, std::ios::binary);
  if(!os.good())
    return false;
  cereal::BinaryInputArchive archive( os );
  archive(poic);
  return true;
}

void SavePOICollection(const POICollection &poic, std::string filename){
  std::ofstream os(filename, std::ios::binary);
  cereal::BinaryOutputArchive archive( os );
  archive(poic);
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
    SavePOICollection(poic, "ztest_poic");
  }

  {
    POICollection poic;
    CHECK(LoadPOICollection(poic,"asdfasfjkwefjewifj")==false);
    CHECK(LoadPOICollection(poic,"ztest_poic")==true);
    CHECK(poic[0].pole.x==a.pole.x);
    CHECK(poic[1].pole.x==b.pole.x);
    CHECK(poic[2].pole.x==c.pole.x);
    CHECK(poic[3].pole.x==d.pole.x);
  }
}




TEST_CASE("POIindex"){
  POICollection poic;
  poic.emplace_back(std::bitset<12>(), Point2D(-93,45).toRadians(), 0);

  SUBCASE("No result"){
    POIindex poii(poic);
    auto result = poii.query(0);
    CHECK(result.size()==0);
  }

  SUBCASE("One result"){
    poic.emplace_back(std::bitset<12>(), Point2D(-93,45).toRadians(), 72.0*DEG_TO_RAD);
    POIindex poii(poic);
    auto result = poii.query(0);
    CHECK(result[0]==1);
  }

  SUBCASE("Check more"){
    poic.emplace_back(std::bitset<12>(), Point2D(-93.1,45.1).toRadians(), 0*DEG_TO_RAD);
    poic.emplace_back(std::bitset<12>(), Point2D(-93.2,45.1).toRadians(), 0*DEG_TO_RAD);
    poic.emplace_back(std::bitset<12>(), Point2D(-93.1,45.2).toRadians(), 0*DEG_TO_RAD);
    poic.emplace_back(std::bitset<12>(), Point2D(-93.2,45.2).toRadians(), 0*DEG_TO_RAD);
    poic.emplace_back(std::bitset<12>(), Point2D(-93.2,45.2).toRadians(), 36*DEG_TO_RAD);
    poic.emplace_back(std::bitset<12>(), Point2D(23,-23.2).toRadians(), 36*DEG_TO_RAD);
    POIindex poii(poic);
    auto result = poii.query(0);
    CHECK(result.size()==4);
    CHECK(result[0]==1);
  }
}
