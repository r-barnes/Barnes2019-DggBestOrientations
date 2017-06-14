#ifndef _land_index_hpp_
#define _land_index_hpp_

#include "Polygon.hpp"
#include "SpIndex.hpp"
#include <string>

class IndexedShapefile {
 public:
  Polygons polys;
  SpIndex sp;
  IndexedShapefile(const std::string shapefile, const std::string layer);
};

#endif