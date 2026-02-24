#ifndef SHORTESTPATHS_H
#define SHORTESTPATHS_H

#include <vector>
#include <unordered_map>

namespace mol2db2 {

/*
 * Single-source shortest paths using Dijkstra's algorithm with priority queue.
 * Calculates distances from a set of initial nodes to all other nodes.
 * 
 * Parameters:
 *   nodes: list of all nodes
 *   edges: dict from nodes to [(neighbor, distance), ...]
 *   startDist: initial distance for starting nodes
 *   initialNodes: list of starting nodes
 * 
 * Returns:
 *   map from node to distance
 */
std::unordered_map<int, int> shortestPaths(
    const std::vector<int>& nodes,
    const std::unordered_map<int, std::vector<std::pair<int, int>>>& edges,
    int startDist,
    const std::vector<int>& initialNodes);

} // namespace mol2db2

#endif // SHORTESTPATHS_H
