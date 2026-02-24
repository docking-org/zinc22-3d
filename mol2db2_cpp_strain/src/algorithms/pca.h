#ifndef PCA_H
#define PCA_H

#include <vector>
#include <array>
#include <utility>

namespace mol2db2 {

/*
 * Principal Component Analysis for 3D points.
 * Uses Eigen for eigenvalue/eigenvector computation.
 */

/* PCA for list of 3D points */
void pca3d(const std::vector<std::array<double, 3>>& pointList,
           std::vector<double>* eigenvalues,
           std::vector<std::vector<double>>* eigenvectors);

/* PCA for N sets of 3D points (list of list of points) */
void pcaN3d(const std::vector<std::vector<std::array<double, 3>>>& pointListList,
            std::vector<double>* eigenvalues,
            std::vector<std::vector<double>>* eigenvectors);

/* Find longest projected direction from PCA */
std::vector<double> findLongestProjectedDirection(
    const std::vector<std::vector<std::array<double, 3>>>& pointListList);

/* Find biggest gap in projected points and split there */
std::pair<std::vector<int>, std::vector<int>> findBiggestGapSplit(
    const std::vector<double>& projectedPts, int overlap = 0);

/* Split based on average (bisective splitting) */
std::pair<std::vector<int>, std::vector<int>> findBisectiveSplit(
    const std::vector<double>& projectedPts, double avgPoint, int overlap = 0);

/* Project points and split into two groups */
std::pair<std::vector<int>, std::vector<int>> findProjectAndSplit(
    const std::vector<std::vector<std::array<double, 3>>>& pointListList,
    bool altSplit = false, int overlap = 0);

/* Find longest direction for 3D points */
std::vector<double> findLongestDirection(
    const std::vector<std::array<double, 3>>& pointList);

/* Sort eigenvectors by eigenvalue (descending) */
std::vector<std::vector<double>> sortDirections(
    const std::vector<std::array<double, 3>>& pointList);

/* Find longest dimension */
double findLongestDimension(const std::vector<std::array<double, 3>>& pointList);

/* Find dimensions in all 3 principal directions */
std::array<double, 3> findDimensions(
    const std::vector<std::array<double, 3>>& pointList);

/* Flatten list of list of points to list of points */
std::vector<std::vector<double>> flatten(
    const std::vector<std::vector<std::array<double, 3>>>& pointListList);

} // namespace mol2db2

#endif // PCA_H
