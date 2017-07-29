#ifndef _poi_hpp_
#define _poi_hpp_

#include "Point.hpp"
#include <limits>
#include <string>
#include <vector>

class Orientation {
 public:
  Orientation()                     = default;
  Orientation(const Orientation &o) = default;
  Orientation(const Point2D &pole0, const double theta0);
  Point2D pole;
  double theta;
};

class OrientationWithStats : public Orientation {
 public:
  std::vector<bool>   overlaps;
  unsigned char       overlap_count;
  double              mindist       = std::numeric_limits<double>::infinity();
  double              maxdist       = -std::numeric_limits<double>::infinity();
  double              avgdist       = 0;
  int                 edge_overlaps = 0;

  OrientationWithStats() = default;
  OrientationWithStats(const Point2D &pole0, const double theta0);
  OrientationWithStats(const Orientation &o);
};

typedef std::vector<Orientation> OCollection;
typedef std::vector<OrientationWithStats> OSCollection;

bool LoadPOICollection(OSCollection &poic, std::string filename);
void SavePOICollection(const OSCollection &poic, std::string filename);

class OrientationGenerator {
 private:
  const double Rearth = 6371;
  long   N;    //Total number of orientations
  long   Nmax; //Number to be generated given the radial_limit.
 public:
  OrientationGenerator(
    const double point_spacingkm, //Approximate distance between points
    const double radial_limit     //Radians from North pole to which orientations should be generated
  );

  Point2D operator()(long i) const;
  long size() const;
};

#endif