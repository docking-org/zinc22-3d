#ifndef DIVISIVE_CLUSTERING_H
#define DIVISIVE_CLUSTERING_H

#include <vector>
#include <array>

namespace mol2db2 {

/*
 * Exception thrown when cluster fails to split (all points project to same point)
 */
class SplitZeroError {
public:
    SplitZeroError(const std::pair<std::vector<int>, std::vector<int>>& splits) 
        : indicesSplit(splits) {}
    
    std::pair<std::vector<int>, std::vector<int>> indicesSplit;
};

/*
 * Helper functions for divisive clustering
 */

/* Get list of points for given indices */
std::vector<std::vector<std::array<double, 3>>> getListForIndices(
    const std::vector<std::vector<std::array<double, 3>>>& pointListList,
    const std::vector<int>& indices);

/* Find longest sublist index */
int findLongestSubList(const std::vector<std::vector<int>>& clusters);

/* Remap split indices onto original list */
std::pair<std::vector<int>, std::vector<int>> findOrigSplitIndices(
    const std::vector<int>& origList,
    const std::pair<std::vector<int>, std::vector<int>>& splitIndices);

/*
 * Main divisive bisective clustering function.
 * Common use is to cluster positions of ligand atoms.
 * Written generally to cluster lists of points.
 * Returns indices into original lists as lists of lists where each sub-list is a cluster.
 * 
 * Parameters:
 *   pointListList: list of list of 3D points to cluster
 *   numClusters: target number of clusters (default 30)
 *   limit: maximum size of largest cluster before stopping (optional)
 *   startClusters: initial cluster assignment (optional)
 *   verbose: print debug info
 *   overlap: allow overlap between clusters
 */
std::vector<std::vector<int>> divisiveClustering(
    const std::vector<std::vector<std::array<double, 3>>>& pointListList,
    int numClusters = 30,
    const int* limit = nullptr,
    const std::vector<std::vector<int>>* startClusters = nullptr,
    bool verbose = false,
    int overlap = 0);

} // namespace mol2db2

#endif // DIVISIVE_CLUSTERING_H
