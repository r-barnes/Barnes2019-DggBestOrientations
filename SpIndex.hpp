#ifndef _SpIndex_hpp_
#define _SpIndex_hpp_

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

class SpIndex {
 private:
  typedef boost::geometry::model::point<int, 2, boost::geometry::cs::cartesian> point;
  typedef boost::geometry::model::box<point> box;
  typedef std::pair<box, int> value;
  typedef boost::geometry::index::rtree< value, boost::geometry::index::rstar<16> > rtree_t;
  std::vector<value> boxes_to_insert;
  rtree_t rtree;

 public:
  SpIndex();
  void addBox(const int xmin, const int ymin, const int xmax, const int ymax, const int id);
  void addBoxDeferred(const int xmin, const int ymin, const int xmax, const int ymax, const int id);
  int queryPoint(const int px, const int py) const;
  void buildIndex();
};

#endif
