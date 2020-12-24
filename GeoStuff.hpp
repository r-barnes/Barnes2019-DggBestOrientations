#ifndef _geostuff_hpp_
#define _geostuff_hpp_

#include "Point.hpp"
#include "Polygon.hpp"
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicLine.hpp>
#include <functional>
#include <memory>

double GeoDistanceFlatEarth(const Point2D &a, const Point2D &b);

double GeoDistanceHaversine(const Point2D &pa, const Point2D &pb);

double GeoDistanceSphere(const Point2D &a, const Point2D &b);
double GeoDistanceEllipsoid(const Point2D &a, const Point2D &b);

Point2D WGS84toEPSG3857(const Point2D &p);

Point3D WGS84toEllipsoidCartesian(const Point2D &p);
Point2D EllipsoidCartesiantoWGS84(const Point3D &p);

Point3D WGS84toSphericalCartesian(const Point2D &p);
Point2D SphericalCartesiantoWGS84(const Point3D &p);

template<class T>
void ToMercator(T &pts);

double EuclideanDistance(const Point3D &a, const Point3D &b);

Polygons ReadShapefile(std::string filename, std::string layername);


///A base class for generating evenly-distributed points between two ends of a
///great- circle or great-ellipse arc. The class is specialized by several types
///of classes which generate arcs on different spheres and ellipsoids.
class GreatCircleGenerator {
 protected:
  double spacing;
  double da;
  unsigned int num_pts;
 public:
  virtual Point2D operator()(const int i) const = 0;
  unsigned int size() const;
  double getSpacing() const;
};

class GeographicLibGreatCircleGenerator : public GreatCircleGenerator  {
 private:
  GeographicLib::Geodesic geod = GeographicLib::Geodesic(10,0);
  GeographicLib::GeodesicLine gline;
 public:
  GeographicLibGreatCircleGenerator(
    double geod_radius,
    double geod_flattening,
    const Point2D &a,
    const Point2D &b,
    const double spacing0
  );
  Point2D operator()(const int i) const override;
};

class EllipsoidalGreatCircleGenerator : public GeographicLibGreatCircleGenerator {
 public:
  EllipsoidalGreatCircleGenerator(const Point2D &a, const Point2D &b, const double spacing0);
};

class SphericalGreatCircleGenerator : public GeographicLibGreatCircleGenerator {
 public:
  SphericalGreatCircleGenerator(const Point2D &a, const Point2D &b, const double spacing0);
};

class SimpleSphericalGreatCircleGenerator : public GreatCircleGenerator {
 private:
  Point3D svec; //Starting vector
  Point3D bvec; //Bearing vector
 public:
  SimpleSphericalGreatCircleGenerator(const Point2D &a, const Point2D &b, const double spacing0);
  Point2D operator()(int i) const override;
};



enum class GreatCircleGeneratorType {
 SPHERICAL,
 ELLIPSOIDAL,
 SIMPLE_SPHERICAL
};

class GreatArcFactory {
 public:
  static std::shared_ptr<GreatCircleGenerator> make(
    const GreatCircleGeneratorType gcgt,
    const Point2D &a,
    const Point2D &b,
    const double spacing0
  );
};



typedef std::function<Point3D(const Point2D&)> TransLLto3D_t;
typedef std::function<Point2D(const Point3D&)> Trans3DtoLL_t;

typedef std::function<double(const Point2D &a, const Point2D &b)> GeoDistance_t;

#endif