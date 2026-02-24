#include "geometry.h"
#include <iostream>
#include <cmath>
#include <cassert>

using namespace mol2db2;

#define ASSERT_NEAR(a, b, tol) assert(fabs((a) - (b)) < (tol))

void test_distL2Squared3() {
    std::array<double, 3> a = {0.0, 0.0, 0.0};
    std::array<double, 3> b = {3.0, 4.0, 0.0};
    
    double dist_sq = distL2Squared3(a, b);
    ASSERT_NEAR(dist_sq, 25.0, 0.0001);
    
    std::cout << "✓ test_distL2Squared3 passed" << std::endl;
}

void test_cross() {
    std::array<double, 3> x = {1.0, 0.0, 0.0};
    std::array<double, 3> y = {0.0, 1.0, 0.0};
    
    std::array<double, 3> result = cross(x, y);
    
    ASSERT_NEAR(result[0], 0.0, 0.0001);
    ASSERT_NEAR(result[1], 0.0, 0.0001);
    ASSERT_NEAR(result[2], 1.0, 0.0001);
    
    std::cout << "✓ test_cross passed" << std::endl;
}

void test_dot() {
    std::vector<double> x = {1.0, 2.0, 3.0};
    std::vector<double> y = {4.0, 5.0, 6.0};
    
    double result = dot(x, y);
    ASSERT_NEAR(result, 32.0, 0.0001);  // 1*4 + 2*5 + 3*6 = 32
    
    std::cout << "✓ test_dot passed" << std::endl;
}

void test_normalizeVector() {
    std::vector<double> v = {3.0, 4.0, 0.0};
    std::vector<double> normalized = normalizeVector(v);
    
    ASSERT_NEAR(normalized[0], 0.6, 0.0001);
    ASSERT_NEAR(normalized[1], 0.8, 0.0001);
    ASSERT_NEAR(normalized[2], 0.0, 0.0001);
    
    // Verify unit length
    double len = sqrt(normalized[0]*normalized[0] + 
                     normalized[1]*normalized[1] + 
                     normalized[2]*normalized[2]);
    ASSERT_NEAR(len, 1.0, 0.0001);
    
    std::cout << "✓ test_normalizeVector passed" << std::endl;
}

void test_getAverage() {
    std::vector<std::array<double, 3>> points = {
        {0.0, 0.0, 0.0},
        {2.0, 0.0, 0.0},
        {0.0, 2.0, 0.0},
        {0.0, 0.0, 2.0}
    };
    
    std::array<double, 3> avg = getAverage(points);
    
    ASSERT_NEAR(avg[0], 0.5, 0.0001);
    ASSERT_NEAR(avg[1], 0.5, 0.0001);
    ASSERT_NEAR(avg[2], 0.5, 0.0001);
    
    std::cout << "✓ test_getAverage passed" << std::endl;
}

void test_rotateAboutLine() {
    // Rotate point (1,0,0) about z-axis by 90 degrees
    std::array<double, 3> a = {0.0, 0.0, 0.0};  // Point on line
    std::array<double, 3> d = {0.0, 0.0, 1.0};  // Direction
    std::array<double, 3> p = {1.0, 0.0, 0.0};  // Point to rotate
    double theta = M_PI / 2.0;  // 90 degrees
    
    std::array<double, 3> result = rotateAboutLine(a, d, p, theta);
    
    ASSERT_NEAR(result[0], 0.0, 0.0001);
    ASSERT_NEAR(result[1], 1.0, 0.0001);
    ASSERT_NEAR(result[2], 0.0, 0.0001);
    
    std::cout << "✓ test_rotateAboutLine passed" << std::endl;
}

void test_withinTolerance() {
    std::array<double, 3> a = {1.0, 2.0, 3.0};
    std::array<double, 3> b = {1.001, 2.001, 3.001};
    std::array<double, 3> c = {1.1, 2.1, 3.1};
    
    assert(withinTolerance(a, b, 0.01));
    assert(!withinTolerance(a, c, 0.01));
    
    std::cout << "✓ test_withinTolerance passed" << std::endl;
}

int main() {
    std::cout << "Running geometry tests..." << std::endl;
    
    test_distL2Squared3();
    test_cross();
    test_dot();
    test_normalizeVector();
    test_getAverage();
    test_rotateAboutLine();
    test_withinTolerance();
    
    std::cout << "\n✓✓✓ All geometry tests passed! ✓✓✓\n" << std::endl;
    return 0;
}
