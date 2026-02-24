#include "shortestpaths.h"
#include "priodict.h"
#include <algorithm>

namespace mol2db2 {

std::unordered_map<int, int> shortestPaths(
    const std::vector<int>& nodes,
    const std::unordered_map<int, std::vector<std::pair<int, int>>>& edges,
    int startDist,
    const std::vector<int>& initialNodes) {
    
    PriorityDictionary<int, int> currentNodes;
    std::unordered_map<int, int> nodeDist;
    
    /* Initialize with starting nodes */
    for (size_t i = 0; i < initialNodes.size(); i++) {
        currentNodes.set(initialNodes[i], startDist);
    }
    
    /* Process nodes in order of distance */
    while (!currentNodes.empty()) {
        int currentNode = currentNodes.smallest();
        int lastDist = currentNodes.get(currentNode);
        currentNodes.pop();
        
        /* Update if we found a shorter path */
        auto it = nodeDist.find(currentNode);
        if (it == nodeDist.end() || lastDist < it->second) {
            nodeDist[currentNode] = lastDist;
            
            /* Update neighbors */
            auto edgeIt = edges.find(currentNode);
            if (edgeIt != edges.end()) {
                const std::vector<std::pair<int, int>>& neighbors = edgeIt->second;
                
                for (size_t i = 0; i < neighbors.size(); i++) {
                    int neighborNode = neighbors[i].first;
                    int nbDist = neighbors[i].second;
                    int newDist = lastDist + nbDist;
                    
                    if (!currentNodes.contains(neighborNode) || 
                        newDist <= currentNodes.get(neighborNode)) {
                        currentNodes.set(neighborNode, newDist);
                    }
                }
            }
        }
    }
    
    return nodeDist;
}

} // namespace mol2db2
