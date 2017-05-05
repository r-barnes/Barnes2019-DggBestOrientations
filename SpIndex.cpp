#include "SpIndex.hpp"
#include <CGAL/Cartesian.h>
#include <CGAL/Segment_tree_k.h>
#include <CGAL/Range_segment_tree_traits.h>

SpIndex::SpIndex(){
  
}

void SpIndex::addBox(
  const int xmin, 
  const int ymin, 
  const int xmax, 
  const int ymax, 
  const int id
){
  intervals_to_add.push_back(Interval(Pure_interval(Key(xmin,ymin), Key(xmax,ymax)),id));
}

void SpIndex::buildIndex(){
  index = Segment_tree_2_type(intervals_to_add.begin(),intervals_to_add.end());
  intervals_to_add.clear();
}

int SpIndex::queryPoint(const int px, const int py) {
  std::list<Interval> output;
  const auto q = Interval(Pure_interval(Key(px,py), Key(px+1,py+1)),0);
  index.window_query(q,std::back_inserter(output));
  if(output.size()==0)
    return -1;
  else
    return output.front().second;
}
