#ifndef _orientation_index_hpp_
#define _orientation_index_hpp_

#include "Orientation.hpp"
#include "nanoflann.hpp"
#include "Solid.hpp"
#include <memory>
#include <limits>
#include <unordered_map>

class OrientationIndex {
 private:
  typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, OrientationIndex > ,
    OrientationIndex,
    3 /* dim */
  > my_kd_tree_t;

  std::unique_ptr<my_kd_tree_t> index;

  std::vector<Point3D> p3ds;
  std::vector<Point2D> p2ds;
  std::vector<unsigned int> pidx;

  void addOrientation(const unsigned int id, const Orientation &o);

  bool vertexInSubdivision(const Point3D &v) const;

 public:
  const double Rearth = 6371;
  static const unsigned int NO_IGNORE = std::numeric_limits<unsigned int>::max();
  
  //mutable std::vector<double> distance_distribution;

  OrientationIndex(const OCollection &orients);
  OrientationIndex(const OSCollection &orients);

  ~OrientationIndex();
  
  void printExtrema() const;

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

  std::vector<std::pair<unsigned int, double> > query(const Point3D &qp, const double distance) const;
  std::unordered_map<unsigned int,double> distancesToNearbyOrientations(const std::vector<Point3D>::const_iterator qvec_begin, const std::vector<Point3D>::const_iterator qvec_end, const unsigned int ignore_pt, const double distance) const;
  std::vector<unsigned int> query(const unsigned int qpn, const double distance) const;
  std::vector<unsigned int> query(const Orientation &o, const double distance) const;

  //void print() const;
};

#endif