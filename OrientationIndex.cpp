#include "OrientationIndex.hpp"
#include "Solid.hpp"
#include "doctest.h"
#include <unordered_map>

OrientationIndex::~OrientationIndex(){
  index->freeIndex();
}

inline size_t OrientationIndex::kdtree_get_point_count() const {
  return p3ds.size();
}



inline bool OrientationIndex::vertexInSubdivision(const Point3D &v) const {
  // return v.z>=0 && v.y>=0;
  (void)v;
  return true;
}


// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
inline double OrientationIndex::kdtree_distance(const double *p1, const size_t idx_p2, size_t /*size*/) const {
  const double d0 = p1[0]-p3ds[idx_p2].x;
  const double d1 = p1[1]-p3ds[idx_p2].y;
  const double d2 = p1[2]-p3ds[idx_p2].z;
  return d0*d0+d1*d1+d2*d2;
}

OrientationIndex::OrientationIndex(
  const OCollection &orients,
  const SolidifyingFunc sf, 
  const int required_ncount0
){  //Cannot be parallelized, otherwise the order of the points might get mixed up
  for(unsigned int oi=0;oi<orients.size();oi++)
    addOrientation(oi, sf(orients[oi]));

  required_ncount = required_ncount0;

  index.reset(new my_kd_tree_t(3 /*dim*/, *this, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) ));
  index->buildIndex();
}

OrientationIndex::OrientationIndex(
  const OSCollection &orients,
  const SolidifyingFunc sf, 
  const int required_ncount0
){  //Cannot be parallelized, otherwise the order of the points might get mixed up
  for(unsigned int oi=0;oi<orients.size();oi++)
    addOrientation(oi, sf(orients[oi]));

  required_ncount = required_ncount0;

  index.reset(new my_kd_tree_t(3 /*dim*/, *this, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) ));
  index->buildIndex();
}

void OrientationIndex::addOrientation(const unsigned int id, const SolidXY &sxy){
  const SolidXYZ sxyz = sxy.toXYZ(Rearth);
  for(unsigned int vi=0;vi<sxyz.v.size();vi++){
    //Choose one quadrant of 3-space. Will have 2-3 members
    if(vertexInSubdivision(sxyz.v[vi])){ 
      p3ds.push_back(sxyz.v[vi]);
      p2ds.push_back(sxy.v[vi]);
      pidx.push_back(id);
    }
  }
}

// Returns the dim'th component of the idx'th point in the class:
// Since this is inlined and the "dim" argument is typically an immediate value, the
//  "if/else's" are actually solved at compile time.
inline double OrientationIndex::kdtree_get_pt(const size_t idx, int dim) const {
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
bool OrientationIndex::kdtree_get_bbox(BBOX& /* bb */) const { return false; }

std::vector<std::pair<unsigned int, double> > OrientationIndex::query(const Point3D &qp, const double distance) const {
  double query_pt[3] = {qp.x,qp.y,qp.z};

  // const size_t num_results = 3000;
  // size_t ret_index[num_results];
  // double out_dist_sqr[num_results];
  // nanoflann::KNNResultSet<double> resultSet(num_results);
  // resultSet.init(ret_index, out_dist_sqr);
  // nanoflann::SearchParams sp{32,0.05,true};
  // index->findNeighbors(resultSet, &query_pt[0], sp);

  // std::vector<std::pair<unsigned int, double> > matches;
  // for(unsigned int i=0;i<resultSet.size();i++)
  //   matches.emplace_back(ret_index[i], out_dist_sqr[i]);

  std::vector<std::pair<size_t, double> > temp;
  nanoflann::SearchParams params;

  index->radiusSearch(query_pt, distance*distance, temp, params);
  std::sort(temp.begin(),temp.end(),[&](const std::pair<size_t, double> &a, const std::pair<size_t, double> &b){return a.second<b.second;});

  std::vector<std::pair<unsigned int, double> > matches;
  for(auto &t: temp)
    matches.emplace_back(t.first, t.second);

  //std::cerr<<"Found "<<matches.size()<<" in radius."<<std::endl;
  //std::cerr<<"Smallest = "<<matches.front().second<<", Largest="<<matches.back().second<<std::endl;

  return matches;
}



std::vector<std::pair<unsigned int,double> > OrientationIndex::queryWithDistance(const unsigned int qpn, const double distance) const {
  //Locate the first point corresponding to qpn, we'll iterate forward to
  //identify the rest
  auto it_pidx = std::lower_bound(pidx.begin(), pidx.end(), qpn);

  std::vector<Point3D>::const_iterator it3d_start = p3ds.begin() + (it_pidx-pidx.begin());
  std::vector<Point3D>::const_iterator it3d_end   = it3d_start;

  //Advance `it_pidx` to the element just past the end of this orientation,
  //advance `it3d_end` along with it.
  while(*it_pidx==qpn){
    it_pidx++;
    it3d_end++;
  }

  const auto ori_dist = distancesToNearbyOrientations(it3d_start,it3d_end,qpn,distance);

  //Put all neighbours which are close enough into the vector
  std::vector<std::pair<unsigned int, double> > ret;
  for(const auto &nd: ori_dist)
    ret.emplace_back(nd.first,nd.second);

  std::sort(ret.begin(), ret.end(), [](const std::pair<unsigned int, double> &a, const std::pair<unsigned int, double> &b){ return a.second<b.second;});

  if(ret.size()>0){
    CHECK(ret.front().second<=ret.back().second);
  }

  return ret;
}

std::vector<unsigned int> OrientationIndex::query(const unsigned int qpn, const double distance) const {
  const auto ori_dist = queryWithDistance(qpn, distance);

  std::vector<unsigned int> closest_n;
  for(const auto &x: ori_dist)
    closest_n.emplace_back(x.first);

  return closest_n;
}



std::vector<std::pair<unsigned int,double> > OrientationIndex::queryWithDistance(const SolidXY &sxy, const double distance) const {
  const SolidXYZ sxyz = sxy.toXYZ(Rearth);
  std::vector<Point3D> qps;
  qps.reserve(20);
  for(const auto v: sxyz.v)
    if(vertexInSubdivision(v))
      qps.push_back(v);

  const auto ori_dist = distancesToNearbyOrientations(qps.begin(), qps.end(), NO_IGNORE, distance);

  //Put all neighbours which are close enough into the vector
  std::vector<std::pair<unsigned int, double> > ret;
  for(const auto &nd: ori_dist)
    ret.emplace_back(nd.first,nd.second);

  std::sort(ret.begin(), ret.end(), [](const std::pair<unsigned int, double> &a, const std::pair<unsigned int, double> &b){ return a.second<b.second;});

  if(ret.size()>0){
    CHECK(ret.front().second<=ret.back().second);
  }

  return ret;
}

std::vector<unsigned int> OrientationIndex::query(const SolidXY &sxy, const double distance) const {
  const auto ori_dist = queryWithDistance(sxy, distance);

  std::vector<unsigned int> closest_n;
  for(const auto &x: ori_dist)
    closest_n.emplace_back(x.first);

  return closest_n;
}



std::unordered_map<unsigned int,double> OrientationIndex::distancesToNearbyOrientations(
  const std::vector<Point3D>::const_iterator qvec_begin,
  const std::vector<Point3D>::const_iterator qvec_end,
  const unsigned int ignore_pt,
  const double distance
) const {
  //For each 3D point of the query POI, find its nearest neighbours in 3-space
  std::vector< std::vector<std::pair<unsigned int, double> > > matches;
  matches.reserve(20); //Can have 2-3 neighbours
  //Iterate through all of the points associated with qpn that we have in the
  //index, finding their nearest neighbours
  for(auto qviter = qvec_begin;qviter!=qvec_end;++qviter)
    matches.push_back(query(*qviter, distance));

  //Look at the neighbours of each point and determine how many times each
  //neighbour is seen total. Ignore qpn itself.
  std::unordered_map<size_t,uint8_t> counts(10000);
  for(const auto &r: matches)
  for(const auto &x: r){
    if(pidx[x.first]!=ignore_pt)
      counts[pidx[x.first]]++;
  }

  //Neighbours seen near two points represent an actual neighbouring
  //configuration, rather than just one that is nearby. Let us remove those
  //configurations which do not fit this criterion.
  for(auto it = counts.begin(); it!=counts.end(); )
    if(it->second!=required_ncount)
      it = counts.erase(it);
    else
      ++it;

  //For those that remain, sum their distances to the query point
  std::unordered_map<unsigned int,double> distances(10000);
  //For each vertex of qpn that falls in the search zone, get the distances from
  //that vertex to all of its nearest neighbours. However, ignore those
  //neighbours which did not also appear near at least one other vertex of qpn.
  for(unsigned int i=0;i<matches.size();i++)
  for(const auto &x: matches[i]){
    if(counts.count(pidx[x.first])==0)
      continue;
    //const auto dist = GeoDistanceHaversine(p2ds[x.first], p2ds[fs+i]);
    const auto dist = x.second;
    if(dist>distance*distance) //Is neighbor too distant? If so, trash the whole orientation (20 km limit)
      distances[pidx[x.first]] += std::numeric_limits<double>::infinity();
    else
      distances[pidx[x.first]] += dist;
  }

  return distances;
}

// void OrientationIndex::print() const {
  // for(unsigned int i=0;i<pidx.size();i++)
    // std::cout<<std::setw(8)<<pidx[i]<<" "<<std::setw(8)<<(p2ds[i].x*RAD_TO_DEG)<<" "<<std::setw(8)<<(p2ds[i].y*RAD_TO_DEG)<<" "<<std::setw(8)<<p3ds[i].x<<" "<<std::setw(8)<<p3ds[i].y<<" "<<std::setw(8)<<p3ds[i].z<<std::endl;
// }



TEST_CASE("OrientationIndex"){
  const double DEG_TO_RAD = M_PI/180.0;
  const double RAD_TO_DEG = 180.0/M_PI;

  OCollection orients;
  orients.emplace_back(Point2D(-93,45).toRadians(), 0);

  SUBCASE("No result"){
    OrientationIndex oidx(orients, OrientationToRegularIcosahedron, 12);
    auto result = oidx.query(0,100);
    CHECK(result.size()==0);
  }

  SUBCASE("One result"){
    orients.emplace_back(Point2D(-93,45).toRadians(), 72.0*DEG_TO_RAD);
    OrientationIndex oidx(orients, OrientationToRegularIcosahedron, 12);
    auto result = oidx.query(0,100);
    CHECK(result.at(0)==1);
  }

  SUBCASE("Partial overlap with no results"){
    orients.emplace_back(Point2D(-93,45).toRadians(), 36.0*DEG_TO_RAD);
    OrientationIndex oidx(orients, OrientationToRegularIcosahedron, 12);
    auto result = oidx.query(0,100);
    CHECK(result.size()==0);
  }

  SUBCASE("Check more"){
    orients.emplace_back(Point2D(-93.1,45.1).toRadians(), 0*DEG_TO_RAD);
    orients.emplace_back(Point2D(-93.2,45.1).toRadians(), 0*DEG_TO_RAD);
    orients.emplace_back(Point2D(-93.1,45.2).toRadians(), 0*DEG_TO_RAD);
    orients.emplace_back(Point2D(-93.2,45.2).toRadians(), 0*DEG_TO_RAD);
    orients.emplace_back(Point2D(-93.0,45.0).toRadians(), 36*DEG_TO_RAD);
    orients.emplace_back(Point2D(23,-23.2).toRadians(), 36*DEG_TO_RAD);
    OrientationIndex oidx(orients, OrientationToRegularIcosahedron, 12);

    {
      auto result = oidx.query(0,100);
      CHECK(result.size()==4);
      CHECK(result.at(0)==1);
    }
    {
      auto result = oidx.query(OrientationToRegularIcosahedron(orients.front()),100);
      CHECK(result.size()==5);
      CHECK(result.at(0)==0);
    }

  }
}


TEST_CASE("OrientationsWithStats"){
  const double DEG_TO_RAD = M_PI/180.0;
  const double RAD_TO_DEG = 180.0/M_PI;

  OSCollection ows;
  ows.emplace_back(Point2D(-93,45).toRadians(), 0);
  ows.emplace_back(Point2D(-93.1,45.1).toRadians(), 0*DEG_TO_RAD);
  ows.emplace_back(Point2D(-93.2,45.1).toRadians(), 0*DEG_TO_RAD);
  ows.emplace_back(Point2D(-93.1,45.2).toRadians(), 0*DEG_TO_RAD);
  ows.emplace_back(Point2D(-93.2,45.2).toRadians(), 0*DEG_TO_RAD);
  ows.emplace_back(Point2D(-93.2,45.2).toRadians(), 36*DEG_TO_RAD);
  ows.emplace_back(Point2D(23,-23.2).toRadians(), 36*DEG_TO_RAD);
  OrientationIndex oidx(ows, OrientationToRegularIcosahedron, 12);
  auto result = oidx.query(0,100);
  CHECK(result.size()==4);
  CHECK(result.at(0)==1);
}
