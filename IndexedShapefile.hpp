#ifndef _land_index_hpp_
#define _land_index_hpp_

#include "Polygon.hpp"
#include "SpIndex.hpp"
#include <vector>
#include <string>

///An indexed shapefile is one which has a quick look-up as to whether a point
///falls within the bounding box of any of the shapefile's polygons. This class
///handles such lookups.
class IndexedShapefile {
 public:
  Polygons          polys;
  std::vector<bool> rectangular; //True if a given polygon is a rectangle
  SpIndex           sp;
  IndexedShapefile(const std::string shapefile, const std::string layer);
  bool isRect(int pid) const;
};

#endif