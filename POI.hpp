#ifndef _poi_hpp_
#define _poi_hpp_

#include "Point.hpp"
#include "nanoflann.hpp"
#include "Icosa.hpp"
#include <memory>
#include <limits>
#include <bitset>
#include <string>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/bitset.hpp>

class POI {
 public:
  static const unsigned int dim = 12;
  std::bitset<dim>    overlaps      = 0;
  Point2D             pole;
  double              rtheta        = 0;
  double              mindist       = std::numeric_limits<double>::infinity();
  double              maxdist       = -std::numeric_limits<double>::infinity();
  double              avgdist       = 0;
  int                 edge_overlaps = 0;
  POI()                             = default;
  POI(const std::bitset<dim> &overlaps0, const Point2D &pole0, double rtheta0);

  template <class Archive>
  void serialize( Archive & ar ){
    ar(
      overlaps,
      pole,
      rtheta,
      mindist,
      maxdist,
      avgdist,
      edge_overlaps
    );
  }
};



class POIindex {
 private:
  typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, POIindex > ,
    POIindex,
    3 /* dim */
  > my_kd_tree_t;

  std::unique_ptr<my_kd_tree_t> index;

  std::vector<Point3D> p3ds;
  std::vector<Point2D> p2ds;
  std::vector<unsigned int> pidx;

 public:
  const double DIST_LIMIT = 2500; //50 km radius

  POIindex(const std::vector<POI> &pois);
  ~POIindex();

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

  std::vector<std::pair<unsigned int, double> > query(const Point3D &qp) const;
  std::vector<unsigned int> query(const unsigned int qpn) const;
};

typedef std::vector<POI> POICollection;

bool LoadPOICollection(POICollection &poic, std::string filename);
void SavePOICollection(const POICollection &poic, std::string filename);

#endif