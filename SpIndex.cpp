#include "SpIndex.hpp"

SpIndex::SpIndex(){
  
}

void SpIndex::addBox(
  const int xmin, 
  const int ymin, 
  const int xmax, 
  const int ymax, 
  const int id
){
  box b(point(xmin,ymin),point(xmax,ymax));
  rtree.insert(std::make_pair(b, id));
}

void SpIndex::addBoxDeferred(
  const int xmin, 
  const int ymin, 
  const int xmax, 
  const int ymax, 
  const int id
){
  box b(point(xmin,ymin),point(xmax,ymax));
  boxes_to_insert.push_back(std::make_pair(b, id));
}

void SpIndex::buildIndex(){
  rtree = rtree_t(boxes_to_insert);
  boxes_to_insert.clear();
  boxes_to_insert.shrink_to_fit();
}

int SpIndex::queryPoint(const int px, const int py) const {
  box query_box(point(px, py), point(px, py));
  std::vector<value> result_s;
  rtree.query(
    boost::geometry::index::intersects(query_box),
    std::back_inserter(result_s)
  );
  if(result_s.size()==0)
    return -1;
  else
    return result_s.front().second;
}
