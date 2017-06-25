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
  return v.z>=0 && v.y>=0;
}


// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
inline double OrientationIndex::kdtree_distance(const double *p1, const size_t idx_p2, size_t /*size*/) const {
  const double d0 = p1[0]-p3ds[idx_p2].x;
  const double d1 = p1[1]-p3ds[idx_p2].y;
  const double d2 = p1[2]-p3ds[idx_p2].z;
  return d0*d0+d1*d1+d2*d2;
}

OrientationIndex::OrientationIndex(const OCollection &orients){  //Cannot be parallelized, otherwise the order of the points might get mixed up
  for(unsigned int oi=0;oi<orients.size();oi++)
    addOrientation(oi, orients[oi]);

  index.reset(new my_kd_tree_t(3 /*dim*/, *this, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) ));
  index->buildIndex();
}

OrientationIndex::OrientationIndex(const OSCollection &orients){  //Cannot be parallelized, otherwise the order of the points might get mixed up
  for(unsigned int oi=0;oi<orients.size();oi++)
    addOrientation(oi, orients[oi]);

  index.reset(new my_kd_tree_t(3 /*dim*/, *this, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) ));
  index->buildIndex();
}

void OrientationIndex::addOrientation(const unsigned int id, const Orientation &o){
  const SolidXY  p2d = SolidXY(o);
  const SolidXYZ p3d = p2d.toXYZ(Rearth);
  for(unsigned int vi=0;vi<p3d.v.size();vi++){
    //Choose one quadrant of 3-space. Will have 2-3 members
    if(vertexInSubdivision(p3d.v[vi])){ 
      p3ds.push_back(p3d.v[vi]);
      p2ds.push_back(p2d.v[vi]);
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

std::vector<std::pair<unsigned int, double> > OrientationIndex::query(const Point3D &qp) const {
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

  index->radiusSearch(query_pt, DIST_LIMIT, temp, params);
  std::sort(temp.begin(),temp.end(),[&](const std::pair<size_t, double> &a, const std::pair<size_t, double> &b){return a.second<b.second;});

  std::vector<std::pair<unsigned int, double> > matches;
  for(auto &t: temp)
    matches.emplace_back(t.first, t.second);

  //std::cerr<<"Found "<<matches.size()<<" in radius."<<std::endl;
  //std::cerr<<"Smallest = "<<matches.front().second<<", Largest="<<matches.back().second<<std::endl;

  return matches;
}



std::vector<unsigned int> OrientationIndex::query(const unsigned int qpn) const {
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

  return query(it3d_start,it3d_end,qpn);
}



std::vector<unsigned int> OrientationIndex::query(const Orientation &o) const {
  const SolidXYZ sxyz = SolidXY(o).toXYZ(Rearth);
  std::vector<Point3D> qps;
  qps.reserve(12);
  for(const auto v: sxyz.v)
    if(vertexInSubdivision(v))
      qps.push_back(v);
  return query(qps.begin(), qps.end(), NO_IGNORE);
}



std::vector<unsigned int> OrientationIndex::query(
  const std::vector<Point3D>::const_iterator qvec_begin,
  const std::vector<Point3D>::const_iterator qvec_end,
  const unsigned int ignore_pt
) const {
  //For each 3D point of the query POI, find its nearest neighbours in 3-space
  std::vector< std::vector<std::pair<unsigned int, double> > > matches;
  matches.reserve(4); //Can have 2-3 neighbours
  //Iterate through all of the points associated with qpn that we have in the
  //index, finding their nearest neighbours
  for(auto qviter = qvec_begin;qviter!=qvec_end;++qviter)
    matches.push_back(query(*qviter));

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
    if(it->second<2)
      it = counts.erase(it);
    else
      ++it;

  //For those that remain, sum their distances to the query point
  std::unordered_map<size_t,double> distances(10000);
  //For each vertex of qpn that falls in the search zone, get the distances from
  //that vertex to all of its nearest neighbours. However, ignore those
  //neighbours which did not also appear near at least one other vertex of qpn.
  for(unsigned int i=0;i<matches.size();i++)
  for(const auto &x: matches[i]){
    if(counts.count(pidx[x.first])==0)
      continue;
    //const auto dist = GeoDistanceHaversine(p2ds[x.first], p2ds[fs+i]);
    const auto dist = x.second;
    if(dist>DIST_LIMIT) //Is neighbor too distant? If so, trash the whole orientation (20 km limit)
      distances[pidx[x.first]] += std::numeric_limits<double>::infinity();
    else
      distances[pidx[x.first]] += dist;
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

  //if(closest_n.size()>30)
  //  closest_n.resize(30);

  return closest_n;
}



TEST_CASE("OrientationIndex"){
  const double DEG_TO_RAD = M_PI/180.0;
  const double RAD_TO_DEG = 180.0/M_PI;

  OCollection orients;
  orients.emplace_back(Point2D(-93,45).toRadians(), 0);

  SUBCASE("No result"){
    OrientationIndex oidx(orients);
    auto result = oidx.query(0);
    CHECK(result.size()==0);
  }

  SUBCASE("One result"){
    orients.emplace_back(Point2D(-93,45).toRadians(), 72.0*DEG_TO_RAD);
    OrientationIndex oidx(orients);
    auto result = oidx.query(0);
    CHECK(result[0]==1);
  }

  SUBCASE("Partial overlap with no results"){
    orients.emplace_back(Point2D(-93,45).toRadians(), 36.0*DEG_TO_RAD);
    OrientationIndex oidx(orients);
    auto result = oidx.query(0);
    CHECK(result.size()==0);
  }

  SUBCASE("Check more"){
    orients.emplace_back(Point2D(-93.1,45.1).toRadians(), 0*DEG_TO_RAD);
    orients.emplace_back(Point2D(-93.2,45.1).toRadians(), 0*DEG_TO_RAD);
    orients.emplace_back(Point2D(-93.1,45.2).toRadians(), 0*DEG_TO_RAD);
    orients.emplace_back(Point2D(-93.2,45.2).toRadians(), 0*DEG_TO_RAD);
    orients.emplace_back(Point2D(-93.2,45.2).toRadians(), 36*DEG_TO_RAD);
    orients.emplace_back(Point2D(23,-23.2).toRadians(), 36*DEG_TO_RAD);
    OrientationIndex oidx(orients);
    auto result = oidx.query(0);
    CHECK(result.size()==4);
    CHECK(result[0]==1);
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
  OrientationIndex oidx(ows);
  auto result = oidx.query(0);
  CHECK(result.size()==4);
  CHECK(result[0]==1);
}
