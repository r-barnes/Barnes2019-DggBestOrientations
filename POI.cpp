#include "POI.hpp"
#include "GeoStuff.hpp"
#include <unordered_map>
#include <iostream>
#include <fstream>

POI::POI(const std::bitset<12> &overlaps, const Point2D &p, double rtheta){
  this->overlaps = overlaps;
  this->pole     = p;
  this->rtheta   = rtheta;
  ico3d          = IcosaXY(p,rtheta).toXYZ();
}

unsigned int POI::size() const {
  return dim;
}




inline size_t POICollection::kdtree_get_point_count() const {
  return POI::dim*pois.size();
}

// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
inline double POICollection::kdtree_distance(const double *p1, const size_t idx_p2, size_t /*size*/) const {
  const double d0 = p1[0]-pois[idx_p2/12].ico3d.v[idx_p2%12].x;
  const double d1 = p1[1]-pois[idx_p2/12].ico3d.v[idx_p2%12].y;
  const double d2 = p1[2]-pois[idx_p2/12].ico3d.v[idx_p2%12].z;
  return d0*d0+d1*d1+d2*d2;
}

// Returns the dim'th component of the idx'th point in the class:
// Since this is inlined and the "dim" argument is typically an immediate value, the
//  "if/else's" are actually solved at compile time.
inline double POICollection::kdtree_get_pt(const size_t idx, int dim) const {
  if (dim==0)
    return pois[idx/12].ico3d.v[idx%12].x;
  else if (dim==1)
    return pois[idx/12].ico3d.v[idx%12].y;
  else
    return pois[idx/12].ico3d.v[idx%12].z;
}

// Optional bounding-box computation: return false to default to a standard bbox computation loop.
//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
template <class BBOX>
bool POICollection::kdtree_get_bbox(BBOX& /* bb */) const { return false; }

void POICollection::buildIndex() {
  if(index!=NULL)
    delete index;
  index = new my_kd_tree_t(3 /*dim*/, *this, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  index->buildIndex();
}

void POICollection::addPOI(const std::bitset<12> &overlaps, const Point2D &p, double rtheta) {
  pois.emplace_back(overlaps, p, rtheta);
}

void POICollection::addPOI(const POI &poi) {
  pois.push_back(poi);
}

POI& POICollection::operator[](unsigned int i){
  return pois[i];
}

const POI& POICollection::operator[](unsigned int i) const{
  return pois[i];
}

unsigned int POICollection::size() const {
  return pois.size();
}

std::vector<size_t> POICollection::query(const Point3D &qp) const {
  double query_pt[3] = {qp.x,qp.y,qp.z};
  const size_t num_results = 10000;
  size_t ret_index[num_results];
  double out_dist_sqr[num_results];
  nanoflann::KNNResultSet<double> resultSet(num_results);
  resultSet.init(ret_index, out_dist_sqr);
  index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));
  std::vector<size_t> temp(ret_index, ret_index+num_results);
  if(size()<num_results)
    temp.resize(size());
  return temp;
}

std::vector<size_t> POICollection::query(const unsigned int qpn) const {
  const auto &qp = pois[qpn];
  const auto qp2d = qp.ico3d.toLatLon();

  //For each 3D point of the query POI, find its nearest neighbours in 3-space
  std::vector< std::vector<size_t> > results(POI::dim);
  for(unsigned int i=0;i<POI::dim;i++)
    results[i] = query(qp.ico3d.v[i]);

  //Look at the neighbours of each point and determine how many times each
  //neighbour is seen total.
  std::unordered_map<size_t,uint8_t> counts;
  for(const auto &r: results)
  for(const auto &x: r){
    if(x/12!=qpn)
      counts[x/12]++;
  }

  //Neighbours seen near all 12 points represent an actual neighbouring
  //configuration, rather than just one that is nearby. Let us remove those
  //configurations which do not fit this criterion.
  for(auto it = counts.begin(); it!=counts.end(); )
    if(it->second!=12)
      it = counts.erase(it);
    else
      ++it;

  //Convert all of the poles to their 2D form so we can get their distances
  std::unordered_map<size_t,IcosaXY> ico2ds;
  for(const auto &x: counts)
    ico2ds[x.first] = pois[x.first].ico3d.toLatLon();

  //For those that remain, sum their distances to the query point
  std::unordered_map<size_t,double> distances;
  for(unsigned int i=0;i<results.size();i++)
  for(const auto &x: results[i]){
    if(counts.count(x/12)==0)
      continue;
    auto dist = GeoDistanceHaversine(ico2ds[x/12].v[x%12], qp2d.v[i]);
    if(dist>200)
      distances[x/12] += std::numeric_limits<double>::infinity();
    else
      distances[x/12] += GeoDistanceHaversine(ico2ds[x/12].v[x%12], qp2d.v[i]);
  }

  //Put all neighbours which are close enough into the vector
  std::vector<size_t> closest_n;
  for(const auto &nd: distances){
    if(nd.second/12<=200)
      closest_n.emplace_back(nd.first);
  }
  std::sort(closest_n.begin(),closest_n.end(), [&](const size_t a, const size_t b){return distances[a]<distances[b];});

  if(closest_n.size()>0)
    std::cerr<<"Closest = "<<distances[closest_n[0]]<<std::endl;

  if(closest_n.size()>10)
    closest_n.resize(10);

  //TODO
  for(auto &x: closest_n)
    std::cerr<<x<<" "<<distances[x]<<"\n";
  std::cerr<<"\n\n"<<std::endl;

  return closest_n;
}

void POICollection::save(std::string filename) const {
  std::ofstream os(filename, std::ios::binary);
  cereal::BinaryOutputArchive archive( os );
  archive(*this);
}

bool POICollection::load(std::string filename) {
  std::ifstream os(filename, std::ios::binary);
  if(!os.good())
    return false;
  cereal::BinaryInputArchive archive( os );
  archive(*this);
  return true;
}