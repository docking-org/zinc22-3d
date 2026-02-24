#include "geometry.h"
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>

namespace mol2db2 {

/* Distance L2 (Euclidean) */
double distL2(const std::vector<double>& x, const std::vector<double>& y) {
    Eigen::Map<const Eigen::VectorXd> ex(x.data(), x.size());
    Eigen::Map<const Eigen::VectorXd> ey(y.data(), y.size());
    return (ex - ey).norm();
}

/* Distance L2 squared */
double distL2Squared(const std::vector<double>& x, const std::vector<double>& y) {
    double sum = 0.0;
    for (size_t i = 0; i < x.size(); i++) {
        double diff = y[i] - x[i];
        sum += diff * diff;
    }
    return sum;
}

/* Distance L2 squared for 3D points */
double distL2Squared3(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    double dx = b[0] - a[0];
    double dy = b[1] - a[1];
    double dz = b[2] - a[2];
    return dx*dx + dy*dy + dz*dz;
}

/* General distance function */
double dist(const std::vector<double>& a, const std::vector<double>& b, 
            const char* metric) {
    if (strcmp(metric, "L2") == 0) {
        double sum = 0.0;
        for (size_t i = 0; i < a.size(); i++) {
            double diff = b[i] - a[i];
            sum += diff * diff;
        }
        return sqrt(sum);
    } else if (strcmp(metric, "LINF") == 0) {
        double maxVal = 0.0;
        for (size_t i = 0; i < a.size(); i++) {
            double newVal = fabs(b[i] - a[i]);
            if (newVal > maxVal) {
                maxVal = newVal;
            }
        }
        return maxVal;
    } else if (strcmp(metric, "L2SQUARED") == 0) {
        double sum = 0.0;
        for (size_t i = 0; i < a.size(); i++) {
            double diff = b[i] - a[i];
            sum += diff * diff;
        }
        return sum;
    } else if (strcmp(metric, "L1") == 0) {
        double sum = 0.0;
        for (size_t i = 0; i < a.size(); i++) {
            sum += fabs(b[i] - a[i]);
        }
        return sum;
    }
    return 0.0;
}

/* Longest and mean distance between points */
void longestAndMeanDist(const std::vector<std::vector<double>>& pts,
                        double* longestDist, double* meanDist) {
    *longestDist = 0.0;
    double sumDists = 0.0;
    int countDists = 0;
    
    for (size_t i = 0; i < pts.size(); i++) {
        for (size_t j = i + 1; j < pts.size(); j++) {
            double thisDist = distL2(pts[i], pts[j]);
            if (thisDist > *longestDist) {
                *longestDist = thisDist;
            }
            sumDists += thisDist;
            countDists++;
        }
    }
    
    *meanDist = (countDists > 0) ? (sumDists / countDists) : 0.0;
}

// add by TEB
void longestAndMeanDist_va3(const std::vector<std::array<double, 3>>& pts,
                        double* longestDist, double* meanDist){

    *longestDist = 0.0;
    double sumDists = 0.0;
    int countDists = 0;
    
    for (size_t i = 0; i < pts.size(); i++) {
        std::vector<double> pts1;
        for( size_t ii = 0; ii < 3; ii++){
            pts1.push_back(pts[i][ii]);
        }
        for (size_t j = i + 1; j < pts.size(); j++) {
            std::vector<double> pts2;
            for( size_t jj = 0; jj < 3; jj++){
                pts2.push_back(pts[j][jj]);
            }
            double thisDist = distL2(pts1, pts2);
            if (thisDist > *longestDist) {
                *longestDist = thisDist;
            }
            sumDists += thisDist;
            countDists++;
        }
    }
    
    *meanDist = (countDists > 0) ? (sumDists / countDists) : 0.0;
}

/* Get angle between two vectors */
double getAngle(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    double ab = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    double aSquared = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    double bSquared = b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
    
    double cosAngle = ab / (sqrt(aSquared) * sqrt(bSquared));
    /* Clamp to [-1, 1] to avoid numerical issues with acos */
    cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
    
    return acos(cosAngle);
}

/* Calculate triangle area using Heron's formula */
double calcTriArea(const std::array<double, 3>& a, 
                   const std::array<double, 3>& b,
                   const std::array<double, 3>& c) {
    std::array<double, 3> dists;
    dists[0] = sqrt(distL2Squared3(a, b));
    dists[1] = sqrt(distL2Squared3(b, c));
    dists[2] = sqrt(distL2Squared3(a, c));
    
    double s = (dists[0] + dists[1] + dists[2]) * 0.5;
    double triArea = sqrt(s * (s - dists[0]) * (s - dists[1]) * (s - dists[2]));
    return triArea;
}

double calcTriAreaList(const std::vector<std::array<double, 3>>& abc) {
    return calcTriArea(abc[0], abc[1], abc[2]);
}

/* Get vector a - b */
std::vector<double> getVector(const std::vector<double>& a, 
                               const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        result[i] = a[i] - b[i];
    }
    return result;
}

/* Get normalized vector a - b */
std::vector<double> getNormalVector(const std::vector<double>& a,
                                     const std::vector<double>& b) {
    return normalizeVector(getVector(a, b));
}

/* Normalize a vector */
std::vector<double> normalizeVector(const std::vector<double>& vector) {
    double total = 0.0;
    for (size_t i = 0; i < vector.size(); i++) {
        total += vector[i] * vector[i];
    }
    total = sqrt(total);
    
    std::vector<double> newVect(vector.size());
    for (size_t i = 0; i < vector.size(); i++) {
        newVect[i] = vector[i] / total;
    }
    return newVect;
}

/* Vector length */
double length(const std::vector<double>& vector) {
    double total = 0.0;
    for (size_t i = 0; i < vector.size(); i++) {
        total += vector[i] * vector[i];
    }
    return sqrt(total);
}

/* Dot product */
double dot(const std::vector<double>& x, const std::vector<double>& y) {
    double result = 0.0;
    for (size_t i = 0; i < x.size(); i++) {
        result += x[i] * y[i];
    }
    return result;
}

/* Cross product for 3D vectors */
std::array<double, 3> cross(const std::array<double, 3>& x,
                             const std::array<double, 3>& y) {
    std::array<double, 3> result;
    result[0] = x[1] * y[2] - x[2] * y[1];
    result[1] = x[2] * y[0] - x[0] * y[2];
    result[2] = x[0] * y[1] - x[1] * y[0];
    return result;
}

/* Get dihedral angle from 4 points */
double getDihedral(const std::array<double, 3>& a,
                   const std::array<double, 3>& b,
                   const std::array<double, 3>& c,
                   const std::array<double, 3>& d) {
    /* Create normal vectors */
    std::vector<double> ab(3), bc(3), cd(3);
    for (int i = 0; i < 3; i++) {
        ab[i] = a[i] - b[i];
        bc[i] = b[i] - c[i];
        cd[i] = c[i] - d[i];
    }
    
    std::vector<double> nab = normalizeVector(ab);
    std::vector<double> nbc = normalizeVector(bc);
    std::vector<double> ncd = normalizeVector(cd);
    
    /* Cross products */
    std::array<double, 3> nab_arr = {nab[0], nab[1], nab[2]};
    std::array<double, 3> nbc_arr = {nbc[0], nbc[1], nbc[2]};
    std::array<double, 3> ncd_arr = {ncd[0], ncd[1], ncd[2]};
    
    std::array<double, 3> cross1_tmp = cross(nab_arr, nbc_arr);
    std::array<double, 3> cross2_tmp = cross(nbc_arr, ncd_arr);
    
    std::vector<double> cross1_vec(3), cross2_vec(3);
    for (int i = 0; i < 3; i++) {
        cross1_vec[i] = cross1_tmp[i];
        cross2_vec[i] = cross2_tmp[i];
    }
    
    std::vector<double> cross1 = normalizeVector(cross1_vec);
    std::vector<double> cross2 = normalizeVector(cross2_vec);
    
    double dotProduct = dot(cross1, cross2);
    dotProduct = std::max(-1.0, std::min(1.0, dotProduct));
    
    double dihedral1 = acos(dotProduct);
    
    /* Figure out +/- direction */
    std::array<double, 3> cross1_arr = {cross1[0], cross1[1], cross1[2]};
    double planeD = calculatePlaneD(cross1_arr, b);
    std::array<double, 4> planeFull = {cross1[0], cross1[1], cross1[2], planeD};
    
    if (!checkPlaneSide(planeFull, d)) {
        dihedral1 = -dihedral1;
    }
    
    return dihedral1;
}

double getDihedralUnited(const std::vector<std::array<double, 3>>& all) {
    return getDihedral(all[0], all[1], all[2], all[3]);
}

/* Rotate point xyz about line d-a by angle theta */
std::array<double, 3> rotateAboutLine(const std::array<double, 3>& aIn,
                                       const std::array<double, 3>& dIn,
                                       const std::array<double, 3>& xyz,
                                       double theta) {
    /* Constrain theta to -pi to +pi */
    while (theta < -M_PI) {
        theta += 2.0 * M_PI;
    }
    while (theta > M_PI) {
        theta -= 2.0 * M_PI;
    }
    
    /* Line direction vector */
    std::array<double, 3> da;
    for (int i = 0; i < 3; i++) {
        da[i] = dIn[i] - aIn[i];
    }
    
    /* Unpack coordinates */
    double a = aIn[0], b = aIn[1], c = aIn[2];
    double u = da[0], v = da[1], w = da[2];
    double x = xyz[0], y = xyz[1], z = xyz[2];
    
    /* Shortcuts */
    std::vector<double> da_vec = {da[0], da[1], da[2]};
    double uvw = length(da_vec);
    double uvw2 = uvw * uvw;
    
    /* Rotation formulas */
    double newX = (a * (v*v + w*w) + u * (-b*v - c*w + u*x + v*y + w*z) +
                   (-a * (v*v + w*w) + u * (b*v + c*w - v*y - w*z) +
                    x * (v*v + w*w)) * cos(theta) +
                   (-c*v + b*w - w*y + v*z) * sin(theta) * uvw) / uvw2;
    
    double newY = (b * (u*u + w*w) + v * (-a*u - c*w + u*x + v*y + w*z) +
                   (-b * (u*u + w*w) + v * (a*u + c*w - u*x - w*z) +
                    y * (u*u + w*w)) * cos(theta) +
                   (c*u - a*w + w*x - u*z) * sin(theta) * uvw) / uvw2;
    
    double newZ = (c * (v*v + u*u) + w * (-a*u - b*v + u*x + v*y + w*z) +
                   (-c * (v*v + u*u) + w * (a*u + b*v - u*x - v*y) +
                    z * (v*v + u*u)) * cos(theta) +
                   (-b*u + a*v - v*x + u*y) * sin(theta) * uvw) / uvw2;
    
    std::array<double, 3> result = {newX, newY, newZ};
    return result;
}

/* Get triangle normal vector */
std::array<double, 3> getTriNormal(const std::array<double, 3>& a,
                                    const std::array<double, 3>& b,
                                    const std::array<double, 3>& c,
                                    bool firstTime) {
    /* Get vectors a-b and c-b */
    std::array<double, 3> vecAB, vecCB;
    for (int i = 0; i < 3; i++) {
        vecAB[i] = a[i] - b[i];
        vecCB[i] = c[i] - b[i];
    }
    
    /* Cross product */
    std::array<double, 3> normal = cross(vecAB, vecCB);
    
    /* Check if all zeros */
    if (firstTime && normal[0] == 0.0 && normal[1] == 0.0 && normal[2] == 0.0) {
        /* Try permuting */
        std::array<double, 3> newNor = getTriNormal(b, c, a, false);
        if (newNor[0] == 0.0 && newNor[1] == 0.0 && newNor[2] == 0.0) {
            std::array<double, 3> lastNo = getTriNormal(c, a, b, false);
            if (lastNo[0] == 0.0 && lastNo[1] == 0.0 && lastNo[2] == 0.0) {
                return lastNo;
            } else {
                std::vector<double> lastNo_vec = {lastNo[0], lastNo[1], lastNo[2]};
                std::vector<double> normalized = normalizeVector(lastNo_vec);
                return {normalized[0], normalized[1], normalized[2]};
            }
        } else {
            std::vector<double> newNor_vec = {newNor[0], newNor[1], newNor[2]};
            std::vector<double> normalized = normalizeVector(newNor_vec);
            return {normalized[0], normalized[1], normalized[2]};
        }
    } else if (!firstTime) {
        return normal;
    } else {
        std::vector<double> normal_vec = {normal[0], normal[1], normal[2]};
        std::vector<double> normalized = normalizeVector(normal_vec);
        return {normalized[0], normalized[1], normalized[2]};
    }
}

std::array<double, 3> getTriNormalList(const std::vector<std::array<double, 3>>& united) {
    return getTriNormal(united[0], united[1], united[2]);
}

/* Average of 3D points */
std::array<double, 3> getAverage(const std::vector<std::array<double, 3>>& listPoints) {
    std::array<double, 3> average = {0.0, 0.0, 0.0};
    for (size_t i = 0; i < listPoints.size(); i++) {
        for (int j = 0; j < 3; j++) {
            average[j] += listPoints[i][j];
        }
    }
    for (int j = 0; j < 3; j++) {
        average[j] /= listPoints.size();
    }
    return average;
}

/* Average of 1D points */
double getAverage1(const std::vector<double>& listPoints) {
    double average = 0.0;
    for (size_t i = 0; i < listPoints.size(); i++) {
        average += listPoints[i];
    }
    average /= listPoints.size();
    return average;
}

/* Average of arbitrary dimension points */
std::vector<double> getAverageArbitraryDimension(
    const std::vector<std::vector<double>>& listPoints, int dimension) {
    std::vector<double> average(dimension, 0.0);
    for (size_t i = 0; i < listPoints.size(); i++) {
        for (int j = 0; j < dimension; j++) {
            average[j] += listPoints[i][j];
        }
    }
    for (int j = 0; j < dimension; j++) {
        average[j] /= listPoints.size();
    }
    return average;
}

/* Calculate plane D parameter */
double calculatePlaneD(const std::array<double, 3>& normal,
                       const std::array<double, 3>& pointOnP) {
    return -(normal[0] * pointOnP[0] + normal[1] * pointOnP[1] + 
             normal[2] * pointOnP[2]);
}

/* Check which side of plane a point is on */
bool checkPlaneSide(const std::array<double, 4>& plane,
                    const std::array<double, 3>& point) {
    double sign = plane[0] * point[0] + plane[1] * point[1] + 
                  plane[2] * point[2] + plane[3];
    return (sign >= 0.0);
}

/* Distance from plane to origin */
double planeDistToOrigin(const std::array<double, 4>& normal) {
    double a = normal[0], b = normal[1], c = normal[2], d = normal[3];
    return d / sqrt(a*a + b*b + c*c);
}

/* Check if two points are within tolerance */
bool withinTolerance(const std::array<double, 3>& pointA,
                     const std::array<double, 3>& pointB,
                     double tolerance) {
    if (fabs(pointA[0] - pointB[0]) < tolerance) {
        if (fabs(pointA[1] - pointB[1]) < tolerance) {
            if (fabs(pointA[2] - pointB[2]) < tolerance) {
                return true;
            }
        }
    }
    return false;
}

/* Perturb triangle slightly for numerical stability */
void perturbTriangle(const std::array<double, 3>& p1,
                     const std::array<double, 3>& p2,
                     const std::array<double, 3>& p3,
                     std::array<double, 3>* p1new,
                     std::array<double, 3>* p2new,
                     std::array<double, 3>* p3new) {
    for (int i = 0; i < 3; i++) {
        (*p1new)[i] = p1[i] + 0.0000001;
        (*p2new)[i] = p2[i] - 0.000001;
        (*p3new)[i] = p3[i] + 0.00001;
    }
}

/* Line-plane intersection using determinants */
bool linePlaneIntersectionNumeric(const std::array<double, 3>& p1,
                                   const std::array<double, 3>& p2,
                                   const std::array<double, 3>& p3,
                                   const std::array<double, 3>& p4,
                                   const std::array<double, 3>& p5,
                                   std::array<double, 3>* result) {
    Eigen::Matrix4d top;
    top << 1.0, 1.0, 1.0, 1.0,
           p1[0], p2[0], p3[0], p4[0],
           p1[1], p2[1], p3[1], p4[1],
           p1[2], p2[2], p3[2], p4[2];
    
    Eigen::Matrix4d bottom;
    bottom << 1.0, 1.0, 1.0, 0.0,
              p1[0], p2[0], p3[0], p5[0] - p4[0],
              p1[1], p2[1], p3[1], p5[1] - p4[1],
              p1[2], p2[2], p3[2], p5[2] - p4[2];
    
    double topDet = top.determinant();
    double botDet = bottom.determinant();
    
    if (topDet == 0.0 || botDet == 0.0) {
        return false;
    }
    
    double t = -topDet / botDet;
    (*result)[0] = p4[0] + (p5[0] - p4[0]) * t;
    (*result)[1] = p4[1] + (p5[1] - p4[1]) * t;
    (*result)[2] = p4[2] + (p5[2] - p4[2]) * t;
    
    return true;
}

/* Check if intersection point is inside triangle */
bool intPointInsideTri(const std::array<double, 3>& p1,
                       const std::array<double, 3>& p2,
                       const std::array<double, 3>& p3,
                       const std::array<double, 3>& intPt) {
    std::array<double, 3> v21, v31, vint1;
    std::array<double, 3> v12, v32, vint3;
    std::array<double, 3> v13, v23, vint2;
    
    for (int i = 0; i < 3; i++) {
        v21[i] = p2[i] - p1[i];
        v31[i] = p3[i] - p1[i];
        vint1[i] = intPt[i] - p1[i];
        
        v12[i] = p1[i] - p3[i];
        v32[i] = p3[i] - p3[i];
        vint3[i] = intPt[i] - p3[i];
        
        v13[i] = p1[i] - p2[i];
        v23[i] = p2[i] - p2[i];
        vint2[i] = intPt[i] - p2[i];
    }
    
    double p2p3ang = getAngle(v21, v31);
    if (p2p3ang < getAngle(v21, vint1) || p2p3ang < getAngle(v31, vint1)) {
        return false;
    }
    
    double p1p2ang = getAngle(v12, v32);
    if (p1p2ang < getAngle(v32, vint3) || p1p2ang < getAngle(v12, vint3)) {
        return false;
    }
    
    double p3p1ang = getAngle(v13, v23);
    if (p3p1ang < getAngle(v13, vint2) || p3p1ang < getAngle(v23, vint2)) {
        return false;
    }
    
    return true;
}

/* Line-sphere intersection */
bool lineSphereIntersection(const std::array<double, 3>& minLine,
                            const std::array<double, 3>& maxLine,
                            const std::array<double, 4>& sphere,
                            std::array<double, 3>* intersection1,
                            std::array<double, 3>* intersection2) {
    /* Move sphere and line so line starts at origin */
    std::array<double, 3> newSphere;
    for (int i = 0; i < 3; i++) {
        newSphere[i] = sphere[i] - minLine[i];
    }
    double radius = sphere[3];
    
    /* Line direction */
    std::array<double, 3> dirLine;
    for (int i = 0; i < 3; i++) {
        dirLine[i] = maxLine[i] - minLine[i];
    }
    std::vector<double> dirLine_vec = {dirLine[0], dirLine[1], dirLine[2]};
    std::vector<double> normalized = normalizeVector(dirLine_vec);
    dirLine[0] = normalized[0];
    dirLine[1] = normalized[1];
    dirLine[2] = normalized[2];
    
    /* Calculate intersection */
    double partA = 0.0, partB = 0.0, partC = 0.0;
    for (int i = 0; i < 3; i++) {
        partA += dirLine[i] * newSphere[i];
        partB += dirLine[i] * dirLine[i];
        partC += newSphere[i] * newSphere[i];
    }
    partC -= radius * radius;
    
    double discriminant = partA*partA - partB*partC;
    if (discriminant < 0.0) {
        return false;
    }
    
    double sqrtDisc = sqrt(discriminant);
    double d1 = (partA + sqrtDisc) / partB;
    double d2 = (partA - sqrtDisc) / partB;
    
    /* Construct output points */
    for (int i = 0; i < 3; i++) {
        (*intersection1)[i] = minLine[i] + dirLine[i] * d1;
        (*intersection2)[i] = minLine[i] + dirLine[i] * d2;
    }
    
    return true;
}

/* Find min/max of spheres */
bool findMinsMaxsSpheres(const std::vector<std::array<double, 4>>& spheres,
                         std::array<double, 3>* mins,
                         std::array<double, 3>* maxs) {
    if (spheres.size() == 0) {
        return false;
    }
    
    for (int xyz = 0; xyz < 3; xyz++) {
        (*mins)[xyz] = spheres[0][xyz] - spheres[0][3];
        (*maxs)[xyz] = spheres[0][xyz] + spheres[0][3];
    }
    
    for (size_t i = 1; i < spheres.size(); i++) {
        for (int xyz = 0; xyz < 3; xyz++) {
            double minVal = spheres[i][xyz] - spheres[i][3];
            double maxVal = spheres[i][xyz] + spheres[i][3];
            if (minVal < (*mins)[xyz]) {
                (*mins)[xyz] = minVal;
            }
            if (maxVal > (*maxs)[xyz]) {
                (*maxs)[xyz] = maxVal;
            }
        }
    }
    
    return true;
}

/* Calculate sphericity */
double calculateSphericity(double area, double volume) {
    return (pow(M_PI, 1.0/3.0) * pow(6.0 * volume, 2.0/3.0)) / area;
}

/* Fix normal zeros */
std::array<double, 3> fixNormalZeros(const std::array<double, 3>& vector) {
    double alpha = 0.0000000000000000001;
    
    if (vector[0] == 0.0 && vector[1] == 0.0 && vector[2] == 0.0) {
        return vector;
    } else if (vector[0] == 0.0 || vector[1] == 0.0 || vector[2] == 0.0) {
        std::array<double, 3> newVec = vector;
        if (vector[0] == 0.0) newVec[0] += alpha;
        if (vector[1] == 0.0) newVec[1] += alpha;
        if (vector[2] == 0.0) newVec[2] += alpha;
        
        std::vector<double> newVec_v = {newVec[0], newVec[1], newVec[2]};
        std::vector<double> normalized = normalizeVector(newVec_v);
        return {normalized[0], normalized[1], normalized[2]};
    } else {
        return vector;
    }
}

} // namespace mol2db2
