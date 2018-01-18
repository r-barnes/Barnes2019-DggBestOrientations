#include "PointCloud.hpp"
#include "doctest.h"
#include <cstdlib>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iostream>

inline size_t PointCloud::kdtree_get_point_count() const {
  return pts.size();
}

// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
inline double PointCloud::kdtree_distance(const double *p1, const size_t idx_p2, size_t /*size*/) const {
  const auto  &p2 = pts.at(idx_p2);
  const double d0 = p1[0]-p2.x;
  const double d1 = p1[1]-p2.y;
  const double d2 = p1[2]-p2.z;
  return d0*d0+d1*d1+d2*d2;
}

// Returns the dim'th component of the idx'th point in the class:
// Since this is inlined and the "dim" argument is typically an immediate value, the
//  "if/else's" are actually solved at compile time.
inline double PointCloud::kdtree_get_pt(const size_t idx, int dim) const {
  if (dim==0)
    return pts.at(idx).x;
  else if (dim==1)
    return pts.at(idx).y;
  else
    return pts.at(idx).z;
}

void PointCloud::newIndex() {
  index.reset(
    new my_kd_tree_t(3 /*dim*/, *this, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) ),
    [](my_kd_tree_t *t){
      std::cerr<<"Freeing PointCloud index."<<std::endl;
      t->freeIndex(); 
    }
  );
}

// Optional bounding-box computation: return false to default to a standard bbox computation loop.
//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
template <class BBOX>
bool PointCloud::kdtree_get_bbox(BBOX& /* bb */) const { return false; }

void PointCloud::buildIndex() {
  std::cerr<<"Building point cloud index..."<<std::endl;
  newIndex();
  index->buildIndex();
}

void PointCloud::addPoint(const Point3D &xyz) {
  pts.push_back(xyz);
}

const Point3D& PointCloud::queryPoint(const Point3D &xyz) const {
  double query_pt[3] = {xyz.x,xyz.y,xyz.z};
  const size_t num_results = 1;
  size_t ret_index;
  double out_dist_sqr;
  nanoflann::KNNResultSet<double> resultSet(num_results);
  resultSet.init(&ret_index, &out_dist_sqr);
  index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));
  return pts.at(ret_index);
}

std::vector<std::pair<unsigned int, double> > PointCloud::queryByDistance(const Point3D &qp, const double distance) const {
  double query_pt[3] = {qp.x,qp.y,qp.z};

  std::vector< std::pair<size_t, double> > temp;
  nanoflann::SearchParams params;

  index->radiusSearch(query_pt, distance*distance, temp, params);
  std::sort(temp.begin(),temp.end(),[&](const std::pair<size_t, double> &a, const std::pair<size_t, double> &b){return a.second<b.second;});

  std::vector<std::pair<unsigned int, double> > matches;
  for(auto &t: temp)
    matches.emplace_back(t.first, t.second);

  return matches;
}

std::vector<std::pair<unsigned int, double> > PointCloud::queryByAnnulus(const Point3D &qp, const double inner_distance, const double outer_distance) const {
  auto outer = queryByDistance(qp,outer_distance); //Points in outer and inner circles
  auto inner = queryByDistance(qp,inner_distance); //Points in inner circle

  assert(inner_distance<outer_distance);

  //Points are returned in sorted order, so we can find the outermost point of
  //the inner ring
  const auto outermost_of_inner = inner.back().first;

  //Now, find the corresponding point in the outer dataset
  auto iouter = outer.begin();
  for(;iouter!=outer.end();iouter++){
    if(iouter->first==outermost_of_inner)
      break;
  }

  //Increment one more, now the iterator points at the beginning of the annulus
  iouter++;

  outer.erase(outer.begin(),iouter);

  return outer;
}

std::vector< std::pair<int,double> >  PointCloud::kNN(const Point3D &qp, const int num_results) const {
  double query_pt[3] = {qp.x,qp.y,qp.z};
  size_t ret_index[num_results];
  double out_dist_sqr[num_results];
  nanoflann::KNNResultSet<double> resultSet(num_results);
  resultSet.init(ret_index, out_dist_sqr);
  index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

  std::vector< std::pair<int,double> > ret;

  for(int i=0;i<num_results;i++)
    ret.emplace_back(ret_index[i],out_dist_sqr[i]);

  std::sort(ret.begin(),ret.end(),[&](const std::pair<int, double> &a, const std::pair<int, double> &b){return a.second<b.second;});

  return ret;
}




TEST_CASE("PointCloud"){
  PointCloud pc;
  auto a = Point3D(1,0,0);
  pc.addPoint(a);
  pc.addPoint(Point3D(0,1,0));
  pc.addPoint(Point3D(0,0,1));
  pc.buildIndex();
  pc.buildIndex(); //Build twice to test success of rebuild
  auto fp = pc.queryPoint(Point3D(0.8,0.1,-0.05));

  //Should return Point
  CHECK(fp.x==a.x);
  CHECK(fp.y==a.y);
  CHECK(fp.z==0);
}
