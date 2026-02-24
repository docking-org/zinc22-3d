#include "pca.h"
#include "geometry.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>

namespace mol2db2 {

void pca3d(const std::vector<std::array<double, 3>>& pointList,
           std::vector<double>* eigenvalues,
           std::vector<std::vector<double>>* eigenvectors) {
    
    /* Build covariance matrix */
    Eigen::Matrix3d matrix = Eigen::Matrix3d::Zero();
    
    std::array<double, 3> avgPt = getAverage(pointList);
    
    for (size_t i = 0; i < pointList.size(); i++) {
        std::array<double, 3> diffs;
        for (int j = 0; j < 3; j++) {
            diffs[j] = pointList[i][j] - avgPt[j];
        }
        
        /* Build upper diagonal */
        matrix(0, 0) += diffs[0] * diffs[0];
        matrix(0, 1) += diffs[0] * diffs[1];
        matrix(0, 2) += diffs[0] * diffs[2];
        matrix(1, 1) += diffs[1] * diffs[1];
        matrix(1, 2) += diffs[1] * diffs[2];
        matrix(2, 2) += diffs[2] * diffs[2];
    }
    
    /* Make symmetric */
    matrix(1, 0) = matrix(0, 1);
    matrix(2, 0) = matrix(0, 2);
    matrix(2, 1) = matrix(1, 2);
    
    /* Compute eigenvalues and eigenvectors */
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(matrix);
    
    eigenvalues->clear();
    eigenvectors->clear();
    
    for (int i = 0; i < 3; i++) {
        eigenvalues->push_back(solver.eigenvalues()[i]);
        std::vector<double> eigvec(3);
        for (int j = 0; j < 3; j++) {
            eigvec[j] = solver.eigenvectors().col(i)[j];
        }
        eigenvectors->push_back(eigvec);
    }
}

std::vector<std::vector<double>> flatten(
    const std::vector<std::vector<std::array<double, 3>>>& pointListList) {
    
    std::vector<std::vector<double>> flattenedList;
    
    for (size_t i = 0; i < pointListList.size(); i++) {
        std::vector<double> flatTemp;
        for (size_t j = 0; j < pointListList[i].size(); j++) {
            for (int k = 0; k < 3; k++) {
                flatTemp.push_back(pointListList[i][j][k]);
            }
        }
        flattenedList.push_back(flatTemp);
    }
    
    return flattenedList;
}

void pcaN3d(const std::vector<std::vector<std::array<double, 3>>>& pointListList,
            std::vector<double>* eigenvalues,
            std::vector<std::vector<double>>* eigenvectors) {
    
    /* Flatten to matrix form */
    std::vector<std::vector<double>> flattenedList = flatten(pointListList);
    
    if (flattenedList.empty() || flattenedList[0].empty()) {
        return;
    }
    
    int length = flattenedList[0].size();
    int nPoints = flattenedList.size();
    
    /* Create Eigen matrix */
    Eigen::MatrixXd data(nPoints, length);
    for (int i = 0; i < nPoints; i++) {
        for (int j = 0; j < length; j++) {
            data(i, j) = flattenedList[i][j];
        }
    }
    
    /* Center the data */
    Eigen::VectorXd mean = data.colwise().mean();
    Eigen::MatrixXd centered = data.rowwise() - mean.transpose();
    
    /* Compute covariance matrix */
    Eigen::MatrixXd cov = (centered.transpose() * centered) / (nPoints - 1);
    
    /* Compute eigenvalues and eigenvectors */
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(cov);
    
    eigenvalues->clear();
    eigenvectors->clear();
    
    for (int i = 0; i < solver.eigenvalues().size(); i++) {
        double eigval = solver.eigenvalues()[i];
        if (eigval < 0.0 && eigval > -1e-10) {
            eigval = 0.0;  /* Numerical noise */
        }
        eigenvalues->push_back(eigval);
        
        std::vector<double> eigvec(length);
        for (int j = 0; j < length; j++) {
            eigvec[j] = solver.eigenvectors().col(i)[j];
        }
        eigenvectors->push_back(eigvec);
    }
}

std::vector<double> findLongestProjectedDirection(
    const std::vector<std::vector<std::array<double, 3>>>& pointListList) {
    
    std::vector<double> eigenvalues;
    std::vector<std::vector<double>> eigenvectors;
    
    pcaN3d(pointListList, &eigenvalues, &eigenvectors);
    
    /* Find maximum eigenvalue */
    double maxVal = 0.0;
    int maxIndex = 0;
    
    for (size_t i = 0; i < eigenvalues.size(); i++) {
        double val = eigenvalues[i];
        if (val < 0.0) val = -val;  /* Handle complex/negative */
        
        if (val > maxVal) {
            maxVal = val;
            maxIndex = i;
        }
    }
    
    return eigenvectors[maxIndex];
}

std::pair<std::vector<int>, std::vector<int>> findBiggestGapSplit(
    const std::vector<double>& projectedPts, int overlap) {
    
    std::vector<double> sortedPts = projectedPts;
    std::sort(sortedPts.begin(), sortedPts.end());
    
    /* Find biggest gap */
    double biggestGap = 0.0;
    int biggestIndex = -1;
    
    for (size_t i = 0; i < sortedPts.size() - 1; i++) {
        double diff = sortedPts[i + 1] - sortedPts[i];
        if (diff > biggestGap) {
            biggestIndex = i;
            biggestGap = diff;
        }
    }
    
    if (biggestIndex == -1) {
        biggestIndex = 0;
    }
    
    double breakPoint = sortedPts[biggestIndex];
    int lowIdx = std::max(biggestIndex - overlap, 0);
    int highIdx = std::min(biggestIndex + overlap, (int)sortedPts.size() - 1);
    double breakPointLow = sortedPts[lowIdx];
    double breakPointHigh = sortedPts[highIdx];
    
    /* Split points */
    std::vector<int> split0, split1;
    for (size_t i = 0; i < projectedPts.size(); i++) {
        if (projectedPts[i] <= breakPointHigh) {
            split0.push_back(i);
        }
        if (projectedPts[i] >= breakPointLow) {
            split1.push_back(i);
        }
    }
    
    return std::make_pair(split0, split1);
}

std::pair<std::vector<int>, std::vector<int>> findBisectiveSplit(
    const std::vector<double>& projectedPts, double avgPoint, int overlap) {
    
    double overlapDist = 0.0;
    
    if (overlap > 0) {
        std::vector<double> dists;
        for (size_t i = 0; i < projectedPts.size(); i++) {
            dists.push_back(fabs(projectedPts[i] - avgPoint));
        }
        std::sort(dists.begin(), dists.end());
        if (overlap < (int)dists.size()) {
            overlapDist = dists[overlap];
        }
    }
    
    std::vector<int> split0, split1;
    for (size_t i = 0; i < projectedPts.size(); i++) {
        if (projectedPts[i] <= avgPoint + overlapDist) {
            split0.push_back(i);
        }
        if (projectedPts[i] >= avgPoint - overlapDist) {
            split1.push_back(i);
        }
    }
    
    /* If second split is empty, try without equality */
    if (split1.empty()) {
        split0.clear();
        split1.clear();
        for (size_t i = 0; i < projectedPts.size(); i++) {
            if (projectedPts[i] < avgPoint + overlapDist) {
                split0.push_back(i);
            }
            if (projectedPts[i] > avgPoint - overlapDist) {
                split1.push_back(i);
            }
        }
    }
    
    return std::make_pair(split0, split1);
}

std::pair<std::vector<int>, std::vector<int>> findProjectAndSplit(
    const std::vector<std::vector<std::array<double, 3>>>& pointListList,
    bool altSplit, int overlap) {
    
    std::vector<double> maxVecRet = findLongestProjectedDirection(pointListList);
    std::vector<std::vector<double>> flattenedList = flatten(pointListList);
    
    /* Project points onto longest direction */
    std::vector<double> projectedPts;
    for (size_t i = 0; i < flattenedList.size(); i++) {
        std::vector<double> projVec(maxVecRet.size());
        for (size_t j = 0; j < maxVecRet.size(); j++) {
            projVec[j] = maxVecRet[j];
        }
        double projectedPt = dot(projVec, flattenedList[i]);
        projectedPts.push_back(projectedPt);
    }
    
    /* Split based on method */
    if (!altSplit) {
        double avgPoint = getAverage1(projectedPts);
        return findBisectiveSplit(projectedPts, avgPoint, overlap);
    } else {
        return findBiggestGapSplit(projectedPts, overlap);
    }
}

std::vector<double> findLongestDirection(
    const std::vector<std::array<double, 3>>& pointList) {
    
    std::vector<double> eigenvalues;
    std::vector<std::vector<double>> eigenvectors;
    
    pca3d(pointList, &eigenvalues, &eigenvectors);
    
    double maxVal = 0.0;
    int maxIndex = 0;
    
    for (size_t i = 0; i < eigenvalues.size(); i++) {
        double val = fabs(eigenvalues[i]);
        if (val > maxVal) {
            maxVal = val;
            maxIndex = i;
        }
    }
    
    return eigenvectors[maxIndex];
}

std::vector<std::vector<double>> sortDirections(
    const std::vector<std::array<double, 3>>& pointList) {
    
    std::vector<double> eigenvalues;
    std::vector<std::vector<double>> eigenvectors;
    
    pca3d(pointList, &eigenvalues, &eigenvectors);
    
    /* Create pairs and sort */
    std::vector<std::pair<double, std::vector<double>>> eigPairs;
    for (size_t i = 0; i < eigenvalues.size(); i++) {
        eigPairs.push_back(std::make_pair(fabs(eigenvalues[i]), eigenvectors[i]));
    }
    
    std::sort(eigPairs.begin(), eigPairs.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });
    
    std::vector<std::vector<double>> result;
    for (size_t i = 0; i < eigPairs.size(); i++) {
        result.push_back(eigPairs[i].second);
    }
    
    return result;
}

double findLongestDimension(const std::vector<std::array<double, 3>>& pointList) {
    std::vector<double> direction = findLongestDirection(pointList);
    
    if (pointList.empty()) {
        return 0.0;
    }
    
    std::vector<double> first(3);
    for (int i = 0; i < 3; i++) {
        first[i] = pointList[0][i];
    }
    
    double minProj = dot(direction, first);
    double maxProj = minProj;
    
    for (size_t i = 1; i < pointList.size(); i++) {
        std::vector<double> pt(3);
        for (int j = 0; j < 3; j++) {
            pt[j] = pointList[i][j];
        }
        
        double proj = dot(direction, pt);
        if (proj < minProj) minProj = proj;
        if (proj > maxProj) maxProj = proj;
    }
    
    return maxProj - minProj;
}

std::array<double, 3> findDimensions(
    const std::vector<std::array<double, 3>>& pointList) {
    
    std::vector<std::vector<double>> directions = sortDirections(pointList);
    std::array<double, 3> dimensions;
    
    for (int d = 0; d < 3; d++) {
        if (pointList.empty()) {
            dimensions[d] = 0.0;
            continue;
        }
        
        std::vector<double> first(3);
        for (int i = 0; i < 3; i++) {
            first[i] = pointList[0][i];
        }
        
        double minProj = dot(directions[d], first);
        double maxProj = minProj;
        
        for (size_t i = 1; i < pointList.size(); i++) {
            std::vector<double> pt(3);
            for (int j = 0; j < 3; j++) {
                pt[j] = pointList[i][j];
            }
            
            double proj = dot(directions[d], pt);
            if (proj < minProj) minProj = proj;
            if (proj > maxProj) maxProj = proj;
        }
        
        dimensions[d] = maxProj - minProj;
    }
    
    return dimensions;
}

} // namespace mol2db2
