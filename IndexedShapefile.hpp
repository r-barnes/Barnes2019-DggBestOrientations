#ifndef _land_index_hpp_
#define _land_index_hpp_

#include "Polygon.hpp"
#include "SpIndex.hpp"
#include <vector>
#include <string>

class IndexedShapefile {
 public:
  Polygons          polys;
  std::vector<bool> rectangular; //True if a given polygon is a rectangle
  SpIndex           sp;
  IndexedShapefile(const std::string shapefile, const std::string layer);
  bool isRect(int pid) const;
};

#endif