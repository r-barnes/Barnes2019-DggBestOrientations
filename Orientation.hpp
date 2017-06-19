#ifndef _poi_hpp_
#define _poi_hpp_

#include "Point.hpp"
#include <limits>
#include <bitset>
#include <string>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/bitset.hpp>

class Orientation {
 public:
  Orientation()                     = default;
  Orientation(const Orientation &o) = default;
  Orientation(const Point2D &pole0, const double theta0);
  Point2D pole;
  double theta;
  template<class Archive>
  void serialize( Archive & ar ){
    ar(pole,theta);
  }
};

class OrientationWithStats : public Orientation {
 public:
  static const unsigned int dim = 12;
  std::bitset<dim>    overlaps      = 0;
  double              mindist       = std::numeric_limits<double>::infinity();
  double              maxdist       = -std::numeric_limits<double>::infinity();
  double              avgdist       = 0;
  int                 edge_overlaps = 0;

  OrientationWithStats() = default;
  OrientationWithStats(const Point2D &pole0, const double theta0);
  OrientationWithStats(const Orientation &o);

  template <class Archive>
  void serialize( Archive & ar ){
    ar(
      pole,
      theta,
      overlaps,
      mindist,
      maxdist,
      avgdist,
      edge_overlaps
    );
  }
};

typedef std::vector<Orientation> OCollection;
typedef std::vector<OrientationWithStats> OSCollection;

bool LoadPOICollection(OSCollection &poic, std::string filename);
void SavePOICollection(const OSCollection &poic, std::string filename);

#endif