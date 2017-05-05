#ifndef _SpIndex_hpp_
#define _SpIndex_hpp_

#include <CGAL/Cartesian.h>
#include <CGAL/Segment_tree_k.h>
#include <CGAL/Range_segment_tree_traits.h>
#include <list>

class SpIndex {
 private:
  typedef CGAL::Cartesian<int> K;
  typedef CGAL::Segment_tree_map_traits_2<K, char> Traits;
  typedef CGAL::Segment_tree_2<Traits> Segment_tree_2_type;
  typedef Traits::Interval Interval;
  typedef Traits::Pure_interval Pure_interval;
  typedef Traits::Key Key;
  typedef std::list<Interval> IList;
  Segment_tree_2_type index;
  IList intervals_to_add;

 public:
  SpIndex();
  void addBox(const int xmin, const int ymin, const int xmax, const int ymax, const int id);
  void buildIndex();
  int queryPoint(const int px, const int py);
};

#endif
