#include "IndexedShapefile.hpp"
#include "Timer.hpp"
#include "GeoStuff.hpp"

IndexedShapefile::IndexedShapefile(const std::string shapefile, const std::string layer){
  {
    std::cerr<<"Reading shapefile '"<<shapefile<<"'..."<<std::endl;
    Timer tmr;
    polys = ReadShapefile(shapefile, layer);
    std::cerr<<"Read "<<polys.size()<<" polygons."<<std::endl;
    std::cerr<<"Time taken = "<<tmr.elapsed()<<std::endl;
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
