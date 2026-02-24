#ifndef BUCKETS2_H
#define BUCKETS2_H

/*
 * Buckets2 — port of Python buckets2.py
 *
 * For each atom in the molecule, bucket() clusters all input conformations
 * by the spatial position of that atom (within `tolerance` Ångströms).
 * It accumulates results into a shared confClusters map keyed by sorted
 * vectors of conformation indices, building up the information needed by
 * countPositionsPython() in Hierarchy.
 */

#include <vector>
#include <array>
#include <map>
#include <unordered_map>
#include <tuple>
#include <cmath>

namespace mol2db2 {

/* Value stored in confClusters */
struct ConfClusterValue {
    int confNum;                                      /* unique id assigned when first seen */
    std::vector<int> atomIds;                         /* atoms (by mol-atom index) in this cluster */
    std::vector<std::array<double, 3>> xyzPositions; /* representative xyz per atom entry */
};

/* Map: sorted vector of conf indices → cluster value */
using ConfClusterMap = std::map<std::vector<int>, ConfClusterValue>;

/* Hash for 3-tuple int key used by the internal hash map */
struct BucketKey3dHash {
    size_t operator()(const std::array<int, 3>& k) const {
        size_t h = std::hash<int>()(k[0]);
        h ^= std::hash<int>()(k[1]) + 0x9e3779b9u + (h << 6) + (h >> 2);
        h ^= std::hash<int>()(k[2]) + 0x9e3779b9u + (h << 6) + (h >> 2);
        return h;
    }
};

class Buckets2 {
public:
    explicit Buckets2(double tolerance);

    /*
     * Process one atom's positions across all conformations.
     *
     * xyzData  – positions of this atom for each input conf (index == conf index)
     * atomId   – index of this atom in the molecule
     * confNum  – next available unique cluster id (threaded through all calls)
     * confClusters – accumulated map updated in-place
     *
     * Returns {npos, newConfNum} where npos is the number of distinct position
     * clusters found for this atom (matches Python bucket() return value),
     * and newConfNum is the updated global counter.
     */
    std::pair<int,int> bucket(const std::vector<std::array<double, 3>>& xyzData,
                              int atomId,
                              int confNum,
                              ConfClusterMap& confClusters);

private:
    double tolerance;
    double tolerance2;
    double bucketsize;   /* tolerance / sqrt(3) — largest cube inscribed in sphere */
};

} // namespace mol2db2

#endif // BUCKETS2_H
