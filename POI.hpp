#ifndef _poi_hpp_
#define _poi_hpp_

#include "Point.hpp"
#include "nanoflann.hpp"
#include "Icosa.hpp"
#include <limits>
#include <bitset>
#include <string>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/bitset.hpp>

class POI {
 private:
  IcosaXYZ _ico3d;
  IcosaXY  _ico2d;
 public:
  static const unsigned int dim = 12;
  std::bitset<dim>    overlaps      = 0;
  Point2D             pole;
  double              rtheta        = 0;
  double              mindist       = std::numeric_limits<double>::infinity();
  double              maxdist       = -std::numeric_limits<double>::infinity();
  double              avgdist       = 0;
  int                 edge_overlaps = 0;
  uint32_t            mindist_n     = 0;
  uint32_t            maxdist_n     = 0;
  uint32_t            avgdist_n     = 0;
  uint32_t            edge_n        = 0;
  POI()                             = default;
  POI(const std::bitset<dim> &overlaps0, const Point2D &pole0, double rtheta0);
  unsigned int size() const;

  const IcosaXYZ& ico3d() const;
  const IcosaXY&  ico2d() const;

  template <class Archive>
  void serialize( Archive & ar ){
    ar(
      overlaps,
      pole,
      rtheta,
      mindist,
      maxdist,
      avgdist,
      edge_overlaps,
      mindist_n,
      maxdist_n,
      avgdist_n,
      edge_n
    );
    _ico2d = IcosaXY(pole,rtheta);
    _ico3d = _ico2d.toXYZ(6371); //Radius of Earth in km
  }
};



class POICollection {
 private:
  typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, POICollection > ,
    POICollection,
    3 /* dim */
  > my_kd_tree_t;

  my_kd_tree_t *index = NULL;

 public:
  std::vector<POI> pois;

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

  void addPOI(const std::bitset<12> &overlaps, const Point2D &p, double rtheta);
  void addPOI(const POI &poi);

  POI& operator[](unsigned int i);
  const POI& operator[](unsigned int i) const;
  unsigned int size() const;

  std::vector<size_t> query(const Point3D &qp) const;
  std::vector<size_t> query(const unsigned int qpn) const;

  template <class Archive>
  void serialize( Archive & ar ) {
    ar(pois);
  }

  void save(std::string filename) const;
  bool load(std::string filename);
};


#endif