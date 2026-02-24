#ifndef FLOYDWARSHALL_H
#define FLOYDWARSHALL_H

#include <vector>
#include <unordered_map>

namespace mol2db2 {

/*
 * Floyd-Warshall all-pairs shortest paths algorithm.
 * Adapted from CLRS algorithms textbook.
 * O(n^3) complexity.
 * 
 * neighbors format: key -> [(neighbor, distance), (neighbor, distance), ...]
 */
void floydWarshall(
    const std::unordered_map<int, std::vector<std::pair<int, int>>>& neighbors,
    std::vector<std::vector<int>>* distances,
    std::unordered_map<int, int>* orderKeys,
    int infinity = 999999999);

} // namespace mol2db2

#endif // FLOYDWARSHALL_H
