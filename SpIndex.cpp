#include "SpIndex.hpp"
#include "doctest.h"

void AddPolygonToSpIndex(const Polygon &poly, SpIndex &sp, const int id){
  const int xmin = poly.minX();
  const int ymin = poly.minY();
  const int xmax = poly.maxX();
  const int ymax = poly.maxY();
  sp.addBoxDeferred(xmin,ymin,xmax,ymax,id);
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

int SpIndex::queryPoint(const Point2D &xy) const {
  box query_box(point(xy.x, xy.y), point(xy.x, xy.y));
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

TEST_CASE("SpIndex"){
  SpIndex sp;
  int id=0;
  for(double y=0;y<1000;y+=100)
  for(double x=0;x<1000;x+=100)
    sp.addBox(x,y,x+100,y+100,id++);

  CHECK(sp.queryPoint(Point2D(350,350))==33);
  CHECK(sp.queryPoint(Point2D(750,550))==57);

  Polygon p;
  p.exterior.emplace_back(1200,1200);
  p.exterior.emplace_back(1200,1300);
  p.exterior.emplace_back(1300,1300);
  p.exterior.emplace_back(1300,1200);

  AddPolygonToSpIndex(p, sp, 347);
  sp.buildIndex();

  CHECK(sp.queryPoint(Point2D(1250,1270))==347);
  CHECK(sp.queryPoint(Point2D(2750,3379))==-1);

  CHECK(p.containsPoint(Point2D(1250,1270)));
  CHECK(p.containsPoint(Point2D(1243,1222)));
  CHECK(!p.containsPoint(Point2D(1194,1222)));
}