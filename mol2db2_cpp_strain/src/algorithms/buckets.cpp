#include "buckets.h"
#include "geometry.h"
#include <cmath>
#include <algorithm>

namespace mol2db2 {

Bucket3d::Bucket3d(const std::vector<std::array<double, 3>>& pointList,
                   double tolerance) {
    this->tolerance2 = tolerance * tolerance;
    this->pointList = pointList;
    
    /* Extract coordinates by dimension */
    for (size_t i = 0; i < pointList.size(); i++) {
        coords[0].push_back(pointList[i][0]);
        coords[1].push_back(pointList[i][1]);
        coords[2].push_back(pointList[i][2]);
    }
    
    /* Find mins and maxs */
    for (int dim = 0; dim < 3; dim++) {
        if (coords[dim].empty()) {
            mins[dim] = 0;
            maxs[dim] = 0;
        } else {
            double minVal = *std::min_element(coords[dim].begin(), coords[dim].end());
            double maxVal = *std::max_element(coords[dim].begin(), coords[dim].end());
            mins[dim] = (int)floor(minVal - tolerance);
            maxs[dim] = (int)ceil(maxVal + tolerance);
        }
    }
    
    /* Initialize buckets */
    for (int dim = 0; dim < 3; dim++) {
        int count = 1 + maxs[dim] - mins[dim];
        buckets[dim].resize(count);
    }
    
    /* Fill buckets */
    for (int dim = 0; dim < 3; dim++) {
        for (size_t i = 0; i < coords[dim].size(); i++) {
            int bucketOne = (int)floor(coords[dim][i] - tolerance) - mins[dim];
            int bucketTwo = (int)floor(coords[dim][i] + tolerance) - mins[dim];
            
            for (int b = bucketOne; b <= bucketTwo; b++) {
                if (b >= 0 && b < (int)buckets[dim].size()) {
                    buckets[dim][b].insert(i);
                }
            }
        }
    }
    
    /* Build possibly nearby points list */
    for (int xCount = 0; xCount < 1 + maxs[0] - mins[0]; xCount++) {
        for (int yCount = 0; yCount < 1 + maxs[1] - mins[1]; yCount++) {
            for (int zCount = 0; zCount < 1 + maxs[2] - mins[2]; zCount++) {
                /* Intersection of three sets */
                std::set<int> newSet;
                
                if (xCount < (int)buckets[0].size() && 
                    yCount < (int)buckets[1].size() && 
                    zCount < (int)buckets[2].size()) {
                    
                    /* Find intersection */
                    std::set<int> temp;
                    std::set_intersection(
                        buckets[0][xCount].begin(), buckets[0][xCount].end(),
                        buckets[1][yCount].begin(), buckets[1][yCount].end(),
                        std::inserter(temp, temp.begin()));
                    
                    std::set_intersection(
                        temp.begin(), temp.end(),
                        buckets[2][zCount].begin(), buckets[2][zCount].end(),
                        std::inserter(newSet, newSet.begin()));
                }
                
                if (!newSet.empty()) {
                    std::vector<int> newVec(newSet.begin(), newSet.end());
                    possiblyNearbyPoints.push_back(newVec);
                }
            }
        }
    }
}

void Bucket3d::getWithinCluster(UnionFind* clusters) {
    for (size_t bucketIdx = 0; bucketIdx < possiblyNearbyPoints.size(); bucketIdx++) {
        std::vector<int>& bucket = possiblyNearbyPoints[bucketIdx];
        
        std::set<int> indicesLeft;
        for (size_t i = 0; i < bucket.size(); i++) {
            indicesLeft.insert(i);
        }
        
        while (!indicesLeft.empty()) {
            auto it = indicesLeft.begin();
            int oneIndex = *it;
            indicesLeft.erase(it);
            
            int oneXyzIndex = bucket[oneIndex];
            
            std::vector<int> thisCluster;
            if ((int)bucket.size() > bigBucket) {
                thisCluster = clusters->getList(oneXyzIndex);
                if ((int)thisCluster.size() >= (int)bucket.size()) {
                    break;  /* Nothing left to union */
                }
            }
            
            const std::array<double, 3>& oneXyz = pointList[oneXyzIndex];
            
            for (size_t twoIndex = 0; twoIndex < bucket.size(); twoIndex++) {
                int twoXyzIndex = bucket[twoIndex];
                
                if ((int)bucket.size() > bigBucket) {
                    /* Check if already in cluster */
                    bool inCluster = false;
                    for (size_t k = 0; k < thisCluster.size(); k++) {
                        if (thisCluster[k] == twoXyzIndex) {
                            inCluster = true;
                            break;
                        }
                    }
                    if (inCluster) {
                        continue;
                    }
                }
                
                const std::array<double, 3>& twoXyz = pointList[twoXyzIndex];
                if (distL2Squared3(oneXyz, twoXyz) < tolerance2) {
                    clusters->unionSets(oneXyzIndex, twoXyzIndex);
                    indicesLeft.erase(twoIndex);
                }
            }
            
            /* Check if all in one cluster */
            if (bucket.size() == pointList.size()) {
                std::vector<std::vector<int>> clusterList = clusters->toLists();
                if (clusterList.size() == 1) {
                    return;  /* All in one cluster, done */
                }
            }
        }
    }
}

std::set<std::pair<int, int>> Bucket3d::getWithin() {
    std::set<std::pair<int, int>> returnPairs;
    
    for (size_t bucketIdx = 0; bucketIdx < possiblyNearbyPoints.size(); bucketIdx++) {
        const std::vector<int>& bucket = possiblyNearbyPoints[bucketIdx];
        
        for (size_t oneIndex = 0; oneIndex < bucket.size(); oneIndex++) {
            int oneXyzIndex = bucket[oneIndex];
            const std::array<double, 3>& oneXyz = pointList[oneXyzIndex];
            
            for (size_t twoIndex = oneIndex + 1; twoIndex < bucket.size(); twoIndex++) {
                int twoXyzIndex = bucket[twoIndex];
                const std::array<double, 3>& twoXyz = pointList[twoXyzIndex];
                
                if (distL2Squared3(oneXyz, twoXyz) < tolerance2) {
                    int first = oneXyzIndex;
                    int second = twoXyzIndex;
                    if (second < first) {
                        std::swap(first, second);
                    }
                    returnPairs.insert(std::make_pair(first, second));
                }
            }
        }
    }
    
    return returnPairs;
}

} // namespace mol2db2
