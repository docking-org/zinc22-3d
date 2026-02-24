#include "buckets2.h"
#include "geometry.h"
#include <algorithm>
#include <cmath>

namespace mol2db2 {

Buckets2::Buckets2(double tolerance)
    : tolerance(tolerance),
      tolerance2(tolerance * tolerance),
      bucketsize(tolerance / std::sqrt(3.0)) {}

std::pair<int,int> Buckets2::bucket(const std::vector<std::array<double, 3>>& xyzData,
                     int atomId,
                     int confNum,
                     ConfClusterMap& confClusters) {
    int n = static_cast<int>(xyzData.size());

    /* Assign each conformation to a spatial bucket.
     * bucketsize = tolerance/sqrt(3) guarantees that any two points in the
     * *same* bucket are within tolerance of each other (max diagonal of cube). */
    std::unordered_map<std::array<int, 3>, std::vector<int>, BucketKey3dHash> buckets;
    buckets.reserve(n);

    for (int i = 0; i < n; i++) {
        std::array<int, 3> key = {
            static_cast<int>(std::floor(xyzData[i][0] / bucketsize)),
            static_cast<int>(std::floor(xyzData[i][1] / bucketsize)),
            static_cast<int>(std::floor(xyzData[i][2] / bucketsize))
        };
        buckets[key].push_back(i);
    }

    /* Track which confs have already been absorbed into a cluster */
    std::vector<bool> processed(n, false);

    int npos = 0;

    for (auto& [bkey, atoms_in_bucket] : buckets) {
        /* Collect unprocessed confs from this bucket */
        std::vector<int> currCluster;
        for (int conf : atoms_in_bucket) {
            if (!processed[conf]) {
                currCluster.push_back(conf);
                processed[conf] = true;
            }
        }
        if (currCluster.empty()) continue;

        const std::array<double, 3>& ref_xyz = xyzData[currCluster[0]];

        /* Check ±2-bucket neighbourhood in all three dimensions. */
        for (int dx = -2; dx <= 2; dx++) {
        for (int dy = -2; dy <= 2; dy++) {
        for (int dz = -2; dz <= 2; dz++) {
            if (dx == 0 && dy == 0 && dz == 0) continue;
            std::array<int, 3> nkey = {bkey[0] + dx, bkey[1] + dy, bkey[2] + dz};
            auto nit = buckets.find(nkey);
            if (nit == buckets.end()) continue;

            for (int conf : nit->second) {
                if (processed[conf]) continue;
                if (distL2Squared3(ref_xyz, xyzData[conf]) <= tolerance2) {
                    currCluster.push_back(conf);
                    processed[conf] = true;
                }
            }
        }}}

        /* Sort cluster to form canonical key */
        std::sort(currCluster.begin(), currCluster.end());

        auto it = confClusters.find(currCluster);
        if (it == confClusters.end()) {
            /* New unique cluster pattern — assign next confNum */
            confClusters[currCluster] = {confNum, {atomId}, {ref_xyz}};
            confNum++;
        } else {
            /* Existing pattern — append this atom to it */
            it->second.atomIds.push_back(atomId);
            it->second.xyzPositions.push_back(ref_xyz);
        }

        npos++;
    }

    return {npos, confNum};
}

} // namespace mol2db2
