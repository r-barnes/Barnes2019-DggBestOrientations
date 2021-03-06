#ifndef _PointCloud_hpp_
#define _PointCloud_hpp_

#include "nanoflann.hpp"
#include "Point.hpp"
#include <stdexcept>
#include <memory>

class PointCloud {
 private:
  typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, PointCloud > ,
    PointCloud,
    3 /* dim */
  > my_kd_tree_t;

  void newIndex();

  std::shared_ptr<my_kd_tree_t> index;

 public:
  std::vector<Point3D> pts;

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const;

  // Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
  inline double kdtree_distance(const double *p1, const size_t idx_p2, size_t /*size*/) const;

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an immediate value, the
  //  "if/else's" are actually solved at compile time.
  inline double kdtree_get_pt(const size_t idx, int dim) const;

  // Optional bounding-box computation: return false to default to a standard bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
  //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /* bb */) const;

  void buildIndex();
  void addPoint(const Point3D &xyz);
  const Point3D& queryPoint(const Point3D &xyz) const;
  std::vector<std::pair<unsigned int, double> > queryByDistance(const Point3D &qp, const double distance) const;
  std::vector<std::pair<unsigned int, double> > queryByAnnulus (const Point3D &qp, const double inner_distance, const double outer_distance) const;
  std::vector< std::pair<int,double> > kNN(const Point3D &qp, const int num_results) const;
};

#endif