#include "floydwarshall.h"
#include <vector>
#include <algorithm>

namespace mol2db2 {

void floydWarshall(
    const std::unordered_map<int, std::vector<std::pair<int, int>>>& neighbors,
    std::vector<std::vector<int>>* distances,
    std::unordered_map<int, int>* orderKeys,
    int infinity) {
    
    int size = neighbors.size();
    
    /* Initialize distance matrix */
    distances->clear();
    distances->resize(size);
    for (int i = 0; i < size; i++) {
        (*distances)[i].resize(size, infinity);
    }
    
    /* Create ordered keys */
    orderKeys->clear();
    std::vector<int> orderedKeys;
    for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
        orderedKeys.push_back(it->first);
    }
    std::sort(orderedKeys.begin(), orderedKeys.end());
    
    /* Map keys to indices */
    for (size_t i = 0; i < orderedKeys.size(); i++) {
        (*orderKeys)[orderedKeys[i]] = i;
    }
    
    /* Set diagonal to 0 */
    for (int i = 0; i < size; i++) {
        (*distances)[i][i] = 0;
    }
    
    /* Initialize from neighbors */
    for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
        int key = it->first;
        const std::vector<std::pair<int, int>>& neighborList = it->second;
        
        for (size_t i = 0; i < neighborList.size(); i++) {
            int neigh = neighborList[i].first;
            int dist = neighborList[i].second;
            (*distances)[(*orderKeys)[key]][(*orderKeys)[neigh]] = dist;
        }
    }
    
    /* Floyd-Warshall algorithm */
    for (int k = 0; k < size; k++) {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                int throughK = (*distances)[i][k] + (*distances)[k][j];
                if (throughK < (*distances)[i][j]) {
                    (*distances)[i][j] = throughK;
                }
            }
        }
    }
}

} // namespace mol2db2
