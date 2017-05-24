#include "PointCloud.hpp"

inline size_t PointCloud::kdtree_get_point_count() const {
  return pts.size();
}

// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
inline double PointCloud::kdtree_distance(const double *p1, const size_t idx_p2, size_t /*size*/) const {
  const double d0 = p1[0]-pts[idx_p2].x;
  const double d1 = p1[1]-pts[idx_p2].y;
  const double d2 = p1[2]-pts[idx_p2].z;
  return d0*d0+d1*d1+d2*d2;
}

// Returns the dim'th component of the idx'th point in the class:
// Since this is inlined and the "dim" argument is typically an immediate value, the
//  "if/else's" are actually solved at compile time.
inline double PointCloud::kdtree_get_pt(const size_t idx, int dim) const {
  if (dim==0)
    return pts[idx].x;
  else if (dim==1)
    return pts[idx].y;
  else
    return pts[idx].z;
}

// Optional bounding-box computation: return false to default to a standard bbox computation loop.
//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
template <class BBOX>
bool PointCloud::kdtree_get_bbox(BBOX& /* bb */) const { return false; }

void PointCloud::buildIndex() {
  if(index!=NULL)
    delete index;
  index = new my_kd_tree_t(3 /*dim*/, *this, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
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
  return pts[ret_index];
}
