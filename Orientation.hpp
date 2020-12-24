#ifndef _poi_hpp_
#define _poi_hpp_

#include "Point.hpp"
#include <limits>
#include <string>
#include <vector>

///Holds the orientation of a polyhedron
class Orientation {
 public:
  Orientation()                     = default;
  Orientation(const Orientation &o) = default;
  Orientation(const Point2D &pole0, const double theta0);
  Point2D pole;    ///< Lat-long positioning of a pole of the polyhedron
  double theta;    ///< Rotation of the polyhedron about said pole
};

///Holds the orientation of a polyhedron along with statistics associated with
///that orientation
class OrientationWithStats : public Orientation {
 public:
  std::vector<bool>   overlaps;       ///< True/False depending on which of the polyhedron's vertices overlap land
  unsigned char       overlap_count;  ///< Number of `overlaps` whose value is true
  double              mindist       = std::numeric_limits<double>::infinity();  ///< Minimum distance of any vertex to a coast
  double              maxdist       = -std::numeric_limits<double>::infinity(); ///< Maximum distance of any vertex to a coast
  double              avgdist       = 0;                                        ///< Average distance of all vertices from coast
  int                 edge_overlaps = 0;                                        ///< Number of sample points along the polyhedron's edges which overlap a coast

  OrientationWithStats() = default;
  OrientationWithStats(const Point2D &pole0, const double theta0);
  OrientationWithStats(const Orientation &o);
};

typedef std::vector<Orientation> OCollection;
typedef std::vector<OrientationWithStats> OSCollection;

bool LoadPOICollection(OSCollection &poic, std::string filename);
void SavePOICollection(const OSCollection &poic, std::string filename);



///Generates a series of equally-spaced orientations about the globe
class OrientationGenerator {
 private:
  const double Rearth = 6371;
  long   N;                       ///< Total number of orientations
  long   Nmax;                    ///< Number to be generated given the radial_limit.
 public:
  OrientationGenerator(
    const double point_spacingkm, ///< Approximate distance between points
    const double radial_limit     ///< Radians from North pole to which orientations should be generated
  );

  Point2D operator()(long i) const; ///< Returns the ith orientation
  long size() const;                ///< Returns the number of orientations to be generated
};

#endif