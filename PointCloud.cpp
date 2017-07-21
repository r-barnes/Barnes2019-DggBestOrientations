#include "PointCloud.hpp"
#include "doctest.h"
#include <cstdlib>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iostream>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

inline size_t PointCloud::kdtree_get_point_count() const {
  return pts.size();
}

// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
inline double PointCloud::kdtree_distance(const double *p1, const size_t idx_p2, size_t /*size*/) const {
  const auto  &p2 = pts.at(idx_p2);
  const double d0 = p1[0]-p2.x;
  const double d1 = p1[1]-p2.y;
  const double d2 = p1[2]-p2.z;
  return d0*d0+d1*d1+d2*d2;
}

// Returns the dim'th component of the idx'th point in the class:
// Since this is inlined and the "dim" argument is typically an immediate value, the
//  "if/else's" are actually solved at compile time.
inline double PointCloud::kdtree_get_pt(const size_t idx, int dim) const {
  if (dim==0)
    return pts.at(idx).x;
  else if (dim==1)
    return pts.at(idx).y;
  else
    return pts.at(idx).z;
}

void PointCloud::newIndex() {
  index.reset(
    new my_kd_tree_t(3 /*dim*/, *this, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) ),
    [](my_kd_tree_t *t){
      std::cerr<<"Freeing PointCloud index."<<std::endl;
      t->freeIndex(); 
    }
  );
}

// Optional bounding-box computation: return false to default to a standard bbox computation loop.
//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
template <class BBOX>
bool PointCloud::kdtree_get_bbox(BBOX& /* bb */) const { return false; }

void PointCloud::buildIndex() {
  newIndex();
  index->buildIndex();
}

void PointCloud::addPoint(const Point3D &xyz) {
  pts.push_back(xyz);
}

const Point3D& PointCloud::queryPoint(const Point3D &xyz) const {
  double query_pt[3] = {xyz.x,xyz.y,xyz.z};
  const size_t num_results = 1;
  size_t ret_index;
  double out_dist_sqr;
  nanoflann::KNNResultSet<double> resultSet(num_results);
  resultSet.init(&ret_index, &out_dist_sqr);
  index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));
  return pts.at(ret_index);
}

void PointCloud::saveToArchive(std::string fileprefix) const {
  {
    std::ofstream os(fileprefix + "-pts.save", std::ios::binary);
    cereal::BinaryOutputArchive archive( os );
    archive(*this);
  }

  std::string fout_cloud_name = fileprefix+"-cloud.save";
  FILE *f = fopen(fout_cloud_name.c_str(),"wb");
  if (!f) throw std::runtime_error("Error writing index file!");
  index->saveIndex(f);
  fclose(f);
}

bool PointCloud::loadFromArchive(std::string fileprefix) {
  std::ifstream os(fileprefix + "-pts.save", std::ios::binary);
  std::string fout_cloud_name = fileprefix+"-cloud.save";
  FILE *f = fopen(fout_cloud_name.c_str(),"rb");
  if(!os.good() || !f)
    return false;

  cereal::BinaryInputArchive archive( os );
  archive(*this);

  newIndex();
  index->loadIndex(f);
  fclose(f);

  return true;
}


TEST_CASE("PointCloud"){
  PointCloud pc;
  auto a = Point3D(1,0,0);
  pc.addPoint(a);
  pc.addPoint(Point3D(0,1,0));
  pc.addPoint(Point3D(0,0,1));
  pc.buildIndex();
  pc.buildIndex(); //Build twice to test success of rebuild
  auto fp = pc.queryPoint(Point3D(0.8,0.1,-0.05));

  //Should return Point
  CHECK(fp.x==a.x);
  CHECK(fp.y==a.y);
  CHECK(fp.z==0);
}

TEST_CASE("PointCloud Save/Load"){
  Point3D qp(rand(),rand(),rand());
  Point3D closest;
  {
    PointCloud pc;
    for(int i=0;i<200;i++)
      pc.addPoint(Point3D(rand(),rand(),rand()));
    pc.buildIndex();
    pc.saveToArchive("test/test_pc_save");
    closest = pc.queryPoint(qp);
  }

  {
    PointCloud pc;
    auto ret = pc.loadFromArchive("asdfasfdjkjdsf");
    CHECK(!ret);
  }

  {
    PointCloud pc;
    pc.loadFromArchive("test/test_pc_save");
    Point3D lclosest = pc.queryPoint(qp);
    CHECK(lclosest.x==closest.x);
    CHECK(lclosest.y==closest.y);
    CHECK(lclosest.z==closest.z);
  }

}