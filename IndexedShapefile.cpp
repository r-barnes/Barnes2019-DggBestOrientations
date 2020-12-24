#include "IndexedShapefile.hpp"
#include "Timer.hpp"
#include "GeoStuff.hpp"

#include <iostream>

///Read a shapefile and build an (R*-tree) index of its polygons
/// @param shapefile Shapefile to read
/// @param layer     Layer from the shapefile
IndexedShapefile::IndexedShapefile(const std::string shapefile, const std::string layer){
  {
    std::cerr<<"Reading shapefile '"<<shapefile<<"'..."<<std::endl;
    Timer tmr;
    polys = ReadShapefile(shapefile, layer);
    std::cerr<<"Read "<<polys.size()<<" polygons."<<std::endl;
    std::cerr<<"Time taken = "<<tmr.elapsed()<<std::endl;
  }

  {
    int rect_count = 0;
    rectangular.resize(polys.size(),false);
    std::cerr<<"Determining boxality of the polygons"<<std::endl;
    for(unsigned int p=0;p<polys.size();p++){
      //Rectangular polygons have size 5 (four points plus one overlap)
      if(polys[p].exterior.size()!=5)
        continue;

      //TODO: Improved rectangularity check

      //Rectangular polygons' first and last points match, so, if this is not
      //the case here, then this cannot be a rectangle.
      if(polys[p].exterior.front().x!=polys[p].exterior.back().x)
        continue;
      if(polys[p].exterior.front().y!=polys[p].exterior.back().y)
        continue;
      rectangular[p] = true;
      rect_count++;
    }
    std::cerr<<"Rectangular polygons = "<<rect_count<<std::endl;
  }

  {
    std::cerr<<"Building index..."<<std::endl;
    Timer tmr;
    for(int64_t i=0;(unsigned)i<polys.size();i++)
      AddPolygonToSpIndex(polys[i], sp, i);
    sp.buildIndex();
    std::cerr<<"Index built."<<std::endl;
    std::cerr<<"Time taken = "<<tmr.elapsed()<<std::endl;
  }
}

///Indicates whether or not a given polygon is a rectangle, as may often happen
///if pre-subdivided polygons are used.
///@param pid Polygon id of the polygon to be checked
///@returns True if the polygon is a rectangle; otherwise, false
bool IndexedShapefile::isRect(int pid) const {
  return rectangular[pid];
}