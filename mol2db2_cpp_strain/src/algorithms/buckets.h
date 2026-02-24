#ifndef BUCKETS_H
#define BUCKETS_H

#include <vector>
#include <set>
#include <array>
#include "unionfind.h"

namespace mol2db2 {

/*
 * Buckets class for fast 3D point overlap detection.
 * Puts all points into buckets in each dimension.
 * Only compares reasonably nearby points.
 * Assumes tolerance << 1 angstrom.
 */
class Bucket3d {
public:
    /*
     * Constructor takes point list and tolerance.
     * Builds bucket structure for fast searching. O(n) time.
     */
    Bucket3d(const std::vector<std::array<double, 3>>& pointList, 
             double tolerance);
    
    /*
     * Souped-up version that puts nearby points into a UnionFind structure.
     * Does every shortcut possible for speed.
     */
    void getWithinCluster(UnionFind* clusters);
    
    /*
     * Returns pairs of points within tolerance.
     * Only compares within buckets.
     * Slower but doesn't require UnionFind structure.
     */
    std::set<std::pair<int, int>> getWithin();
    
private:
    static const int bigBucket = 100;  /* Arbitrary huge bucket size */
    
    double tolerance2;  /* Squared tolerance */
    std::vector<std::array<double, 3>> pointList;
    std::array<std::vector<double>, 3> coords;  /* Coordinates by dimension */
    std::array<int, 3> mins;
    std::array<int, 3> maxs;
    std::array<std::vector<std::set<int>>, 3> buckets;  /* Buckets per dimension */
    std::vector<std::vector<int>> possiblyNearbyPoints;  /* List of point sets */
};

} // namespace mol2db2

#endif // BUCKETS_H
