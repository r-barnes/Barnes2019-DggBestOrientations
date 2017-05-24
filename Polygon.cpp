#include "Polygon.hpp"
#include <cmath>

double Polygon::minX() const {
  double minx=std::numeric_limits<double>::infinity();
  for(const auto &p: exterior)
    minx = std::min(p.x,minx);
  return minx;
}

double Polygon::maxX() const {
  double maxx=-std::numeric_limits<double>::infinity();
  for(const auto &p: exterior)
    maxx = std::max(p.x,maxx);
  return maxx;
}

double Polygon::minY() const {
  double miny=std::numeric_limits<double>::infinity();
  for(const auto &p: exterior)
    miny=std::min(p.y,miny);
  return miny;
}

double Polygon::maxY() const {
  double maxy=-std::numeric_limits<double>::infinity();
  for(const auto &p: exterior)
    maxy=std::max(p.y,maxy);
  return maxy;
}

void Polygon::toRadians() {
  for(auto &p: exterior)
    p.toRadians();
}

void Polygon::toDegrees() {
  for(auto &p: exterior)
    p.toDegrees();
}

//Info: https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
//Info: https://stackoverflow.com/questions/8721406/how-to-determine-if-a-point-is-inside-a-2d-convex-polygon/23223947#23223947
//Info: http://stackoverflow.com/a/2922778/752843
bool Polygon::containsPoint(const Point2D &xy) const {
  unsigned int i, j;
  int c = 0;
  for (i = 0, j = exterior.size()-1; i < exterior.size(); j = i++) {
    if ( ((exterior[i].y>xy.y) != (exterior[j].y>xy.y)) &&
     (xy.x < (exterior[j].x-exterior[i].x) * (xy.y-exterior[i].y) / (exterior[j].y-exterior[i].y) + exterior[i].x) )
       c = !c;
  }
  return c;
}
