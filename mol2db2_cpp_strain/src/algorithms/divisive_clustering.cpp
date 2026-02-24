#include "divisive_clustering.h"
#include "pca.h"
#include <iostream>
#include <stdexcept>

namespace mol2db2 {

std::vector<std::vector<std::array<double, 3>>> getListForIndices(
    const std::vector<std::vector<std::array<double, 3>>>& pointListList,
    const std::vector<int>& indices) {
    
    std::vector<std::vector<std::array<double, 3>>> newListList;
    for (size_t i = 0; i < indices.size(); i++) {
        newListList.push_back(pointListList[indices[i]]);
    }
    return newListList;
}

int findLongestSubList(const std::vector<std::vector<int>>& clusters) {
    int longestIndex = -1;
    
    for (size_t i = 0; i < clusters.size(); i++) {
        if (longestIndex == -1) {
            longestIndex = i;
        } else if (clusters[longestIndex].size() < clusters[i].size()) {
            longestIndex = i;
        }
    }
    
    return longestIndex;
}

std::pair<std::vector<int>, std::vector<int>> findOrigSplitIndices(
    const std::vector<int>& origList,
    const std::pair<std::vector<int>, std::vector<int>>& splitIndices) {
    
    std::vector<int> newSplit0, newSplit1;
    
    for (size_t i = 0; i < splitIndices.first.size(); i++) {
        newSplit0.push_back(origList[splitIndices.first[i]]);
    }
    
    for (size_t i = 0; i < splitIndices.second.size(); i++) {
        newSplit1.push_back(origList[splitIndices.second[i]]);
    }
    
    return std::make_pair(newSplit0, newSplit1);
}

std::vector<std::vector<int>> divisiveClustering(
    const std::vector<std::vector<std::array<double, 3>>>& pointListList,
    int numClusters,
    const int* limit,
    const std::vector<std::vector<int>>* startClusters,
    bool verbose,
    int overlap) {
    
    /* Initialize clusters */
    std::vector<std::vector<int>> clusters;
    
    if (startClusters != nullptr) {
        clusters = *startClusters;
    } else {
        /* Start with all points in one cluster */
        std::vector<int> allIndices;
        for (size_t i = 0; i < pointListList.size(); i++) {
            allIndices.push_back(i);
        }
        clusters.push_back(allIndices);
    }
    
    /* Keep splitting until we have enough clusters */
    while ((int)clusters.size() < numClusters) {
        /* Find largest cluster */
        int biggestClusterIndex = findLongestSubList(clusters);
        
        if (biggestClusterIndex == -1) {
            break;
        }
        
        std::vector<int> biggestCluster = clusters[biggestClusterIndex];
        
        if (biggestCluster.size() == 1) {
            break;  /* All clusters are singletons */
        }
        
        if (limit != nullptr && *limit > (int)biggestCluster.size()) {
            break;  /* All clusters small enough */
        }
        
        /* Get points for this cluster */
        std::vector<std::vector<std::array<double, 3>>> clusterToSplit = 
            getListForIndices(pointListList, biggestCluster);
        
        try {
            /* Try to split */
            std::pair<std::vector<int>, std::vector<int>> splitIndices = 
                findProjectAndSplit(clusterToSplit, false, overlap);
            
            if (splitIndices.first.empty() || splitIndices.second.empty()) {
                throw SplitZeroError(splitIndices);
            }
            
            /* Remap to original indices */
            std::pair<std::vector<int>, std::vector<int>> origSplits = 
                findOrigSplitIndices(biggestCluster, splitIndices);
            
            /* Remove old cluster and add two new ones */
            clusters.erase(clusters.begin() + biggestClusterIndex);
            clusters.push_back(origSplits.first);
            clusters.push_back(origSplits.second);
            
        } catch (const SplitZeroError& e) {
            if (verbose) {
                std::cout << "projection/split problem during clustering. quitting with "
                         << clusters.size() << " clusters, which should be enough for anybody"
                         << std::endl;
            }
            break;
        } catch (const std::exception& e) {
            if (verbose) {
                std::cout << "convergence problem during clustering. quitting with "
                         << clusters.size() << " clusters, which should be enough for anybody"
                         << std::endl;
            }
            break;
        }
        
        if (verbose) {
            std::cout << "size of each cluster" << std::endl;
            for (size_t i = 0; i < clusters.size(); i++) {
                std::cout << clusters[i].size() << " ";
            }
            std::cout << std::endl;
        }
    }
    
    return clusters;
}

} // namespace mol2db2
