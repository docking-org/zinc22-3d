#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <array>
#include <cmath>

namespace mol2db2 {

/* Distance functions */
double distL2(const std::vector<double>& x, const std::vector<double>& y);
double distL2Squared(const std::vector<double>& x, const std::vector<double>& y);
double distL2Squared3(const std::array<double, 3>& a, const std::array<double, 3>& b);

/* General distance function with different metrics */
double dist(const std::vector<double>& a, const std::vector<double>& b, 
            const char* metric);

/* Longest and mean distance between all pairs of points */
void longestAndMeanDist(const std::vector<std::vector<double>>& pts,
                        double* longestDist, double* meanDist);

// added by TEB
void longestAndMeanDist_va3(const std::vector<std::array<double, 3>>& pts,
                        double* longestDist, double* meanDist);

/* Angle calculations */
double getAngle(const std::array<double, 3>& a, const std::array<double, 3>& b);

/* Triangle area calculations */
double calcTriArea(const std::array<double, 3>& a, 
                   const std::array<double, 3>& b,
                   const std::array<double, 3>& c);
double calcTriAreaList(const std::vector<std::array<double, 3>>& abc);

/* Vector operations */
std::vector<double> getVector(const std::vector<double>& a, 
                               const std::vector<double>& b);
std::vector<double> getNormalVector(const std::vector<double>& a,
                                     const std::vector<double>& b);
std::vector<double> normalizeVector(const std::vector<double>& vector);
double length(const std::vector<double>& vector);
double dot(const std::vector<double>& x, const std::vector<double>& y);
std::array<double, 3> cross(const std::array<double, 3>& x,
                             const std::array<double, 3>& y);

/* Dihedral angle calculations */
double getDihedral(const std::array<double, 3>& a,
                   const std::array<double, 3>& b,
                   const std::array<double, 3>& c,
                   const std::array<double, 3>& d);
double getDihedralUnited(const std::vector<std::array<double, 3>>& all);

/* Rotation about a line */
std::array<double, 3> rotateAboutLine(const std::array<double, 3>& aIn,
                                       const std::array<double, 3>& dIn,
                                       const std::array<double, 3>& xyz,
                                       double theta);

/* Triangle normal vector */
std::array<double, 3> getTriNormal(const std::array<double, 3>& a,
                                    const std::array<double, 3>& b,
                                    const std::array<double, 3>& c,
                                    bool firstTime = true);
std::array<double, 3> getTriNormalList(const std::vector<std::array<double, 3>>& united);

/* Averaging functions */
std::array<double, 3> getAverage(const std::vector<std::array<double, 3>>& listPoints);
double getAverage1(const std::vector<double>& listPoints);
std::vector<double> getAverageArbitraryDimension(
    const std::vector<std::vector<double>>& listPoints, int dimension);

/* Plane calculations */
double calculatePlaneD(const std::array<double, 3>& normal,
                       const std::array<double, 3>& pointOnP);
bool checkPlaneSide(const std::array<double, 4>& plane,
                    const std::array<double, 3>& point);
double planeDistToOrigin(const std::array<double, 4>& normal);

/* Tolerance checking */
bool withinTolerance(const std::array<double, 3>& pointA,
                     const std::array<double, 3>& pointB,
                     double tolerance);

/* Triangle perturbation for numerical stability */
void perturbTriangle(const std::array<double, 3>& p1,
                     const std::array<double, 3>& p2,
                     const std::array<double, 3>& p3,
                     std::array<double, 3>* p1new,
                     std::array<double, 3>* p2new,
                     std::array<double, 3>* p3new);

/* Line-plane intersection */
bool linePlaneIntersectionNumeric(const std::array<double, 3>& p1,
                                   const std::array<double, 3>& p2,
                                   const std::array<double, 3>& p3,
                                   const std::array<double, 3>& p4,
                                   const std::array<double, 3>& p5,
                                   std::array<double, 3>* result);

/* Point inside triangle check */
bool intPointInsideTri(const std::array<double, 3>& p1,
                       const std::array<double, 3>& p2,
                       const std::array<double, 3>& p3,
                       const std::array<double, 3>& intPt);

/* Line-sphere intersection */
bool lineSphereIntersection(const std::array<double, 3>& minLine,
                            const std::array<double, 3>& maxLine,
                            const std::array<double, 4>& sphere,
                            std::array<double, 3>* intersection1,
                            std::array<double, 3>* intersection2);

/* Find min/max of spheres */
bool findMinsMaxsSpheres(const std::vector<std::array<double, 4>>& spheres,
                         std::array<double, 3>* mins,
                         std::array<double, 3>* maxs);

/* Sphericity calculation */
double calculateSphericity(double area, double volume);

/* Helper function for fixing normal zeros */
std::array<double, 3> fixNormalZeros(const std::array<double, 3>& vector);

} // namespace mol2db2

#endif // GEOMETRY_H
