// 2D nearestPoint on PLC unit test.

#include <iostream>
#include <vector>
#include <set>
#include <array>
#include <stdlib.h>
#include <limits>
#include <sstream>
#include <ctime>

#include "polytope.hh"
#include "intersect.hh"
#include "polytope_test_utilities.hh"

using namespace std;

//------------------------------------------------------------------------------
// Helper to check if a point is in the result set
//------------------------------------------------------------------------------
bool containsPoint(const vector<double>& result,
                   const double* point,
                   const double tol) {
  for (int i = 0; i < result.size()/2; ++i) {
    double dx = result[2*i] - point[0];
    double dy = result[2*i+1] - point[1];
    if (dx*dx + dy*dy < tol*tol) {
      return true;
    }
  }
  return false;
}

//------------------------------------------------------------------------------
// The test itself.
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

  // Test tolerance
  double tol = 1.E-9;

  // Create a square PLC with one hole.
  // Outer square: [-5, 5]^2
  // Inner square (hole): [-1, 1]^2
  const unsigned numVertices = 8;
  double vertices[16] = {-5.0, -5.0,
                          5.0, -5.0,
                          5.0,  5.0,
                         -5.0,  5.0,
                         -1.0, -1.0,
                          1.0, -1.0,
                          1.0,  1.0,
                         -1.0,  1.0};
  polytope::PLC<2> plc;
  plc.facets.resize(4);
  plc.holes.resize(1);
  plc.holes[0].resize(4);
  for (unsigned i = 0; i != 4; ++i) {
    plc.facets[i].resize(2);
    plc.facets[i][0] = i;
    plc.facets[i][1] = (i + 1) % 4;
    plc.holes[0][i].resize(2);
    plc.holes[0][i][0] = 4 + i;
    plc.holes[0][i][1] = 4 + (i + 1) % 4;
  }
  
  { // Test 1: 1 intersection with outside boundary
    double p1[2] = {-10.0, 0.0};
    double p2[2] = {-4.0, 0.0};
    std::vector<double> result;
    const unsigned nint = intersect(p1, p2, numVertices, vertices, plc, result);
    POLY_CHECK2(result.size() == 1, "Test 1: size=" << result.size());
    double expected[2] = {-5.0, 0.0};
    POLY_CHECK2(containsPoint(result, expected, tol), "Test 1: missing expected intersection");
    cerr << "Test 1 passed: 1 intersection with outside boundary" << endl;
  }

  { // Test 2: 2 intersections - entry and exit through opposite sides
    double p1[2] = {-6.0, 2.0};
    double p2[2] = {6.0, 2.0};
    std::vector<double> result;
    const unsigned nint = intersect(p1, p2, numVertices, vertices, plc, result);
    POLY_CHECK2(result.size() == 2, "Test 2: size=" << result.size());
    double expected1[2] = {-5.0, 2.0};
    double expected2[2] = {5.0, 2.0};
    POLY_CHECK2(containsPoint(result, expected1, tol), "Test 2: missing entry intersection");
    POLY_CHECK2(containsPoint(result, expected2, tol), "Test 2: missing exit intersection");
    cerr << "Test 2 passed: 2 intersections through opposite sides" << endl;
  }

  { // Test 3: Intersection with corner vertex
    double p1[2] = {4.0, -6.0};
    double p2[2] = {6.0, -4.0};
    std::vector<double> result;
    const unsigned nint = intersect(p1, p2, numVertices, vertices, plc, result);
    POLY_CHECK2(result.size() == 1, "Test 3: size=" << result.size());
    double expected[2] = {5.0, -5.0};
    POLY_CHECK2(containsPoint(result, expected, tol), "Test 3: missing corner intersection");
    cerr << "Test 3 passed: intersection at corner" << endl;
  }

  { // Test 4: Ray passing through hole (4 intersections total)
    double p1[2] = {0.0, -6.0};
    double p2[2] = {0.0, 6.0};
    std::vector<double> result;
    const unsigned nint = intersect(p1, p2, numVertices, vertices, plc, result);
    POLY_CHECK2(result.size() == 4, "Test 4: size=" << result.size() << " (expected 4)");
    double expected[4][2] = {
      {0.0, -5.0},  // outer bottom
      {0.0, -1.0},  // hole bottom
      {0.0, 1.0},   // hole top
      {0.0, 5.0}    // outer top
    };
    for (int i = 0; i < 4; ++i) {
      POLY_CHECK2(containsPoint(result, expected[i], tol),
                  "Test 4: missing intersection " << i);
    }
    cerr << "Test 4 passed: 4 intersections through hole" << endl;
  }

  { // Test 5: Ray along edge
    double p1[2] = {5.0, -6.0};
    double p2[2] = {5.0, 6.0};
    std::vector<double> result;
    const unsigned nint = intersect(p1, p2, numVertices, vertices, plc, result);
    POLY_CHECK2(result.size() == 2, "Test 5: size=" << result.size() << " (expected 2)");
    double expected1[2] = {5.0, -5.0};
    double expected2[2] = {5.0, 5.0};
    POLY_CHECK2(containsPoint(result, expected1, tol), "Test 5: missing bottom intersection");
    POLY_CHECK2(containsPoint(result, expected2, tol), "Test 5: missing top intersection");
    cerr << "Test 5 passed: 2 intersections along edge" << endl;
  }

  { // Test 6: Diagonal ray through square (should hit 4 intersections with hole)
    double p1[2] = {-6.0, -6.0};
    double p2[2] = {6.0, 6.0};
    std::vector<double> result;
    const unsigned nint = intersect(p1, p2, numVertices, vertices, plc, result);
    POLY_CHECK2(result.size() == 4, "Test 6: size=" << result.size() << " (expected 4)");
    double expected[4][2] = {
      {-5.0, -5.0},  // outer entry
      {-1.0, -1.0},  // hole entry
      {1.0, 1.0},    // hole exit
      {5.0, 5.0}     // outer exit
    };
    for (int i = 0; i < 4; ++i) {
      POLY_CHECK2(containsPoint(result, expected[i], tol),
                  "Test 6: missing diagonal intersection " << i);
    }
    cerr << "Test 6 passed: 4 diagonal intersections" << endl;
  }

  { // Test 7: Ray parallel to edge (no intersection)
    double p1[2] = {-6.0, 0.0};
    double p2[2] = {-6.0, 3.0};
    std::vector<double> result;
    const unsigned nint = intersect(p1, p2, numVertices, vertices, plc, result);
    POLY_CHECK2(result.size() == 0, "Test 7: size=" << result.size() << " (expected 0)");
    cerr << "Test 7 passed: no intersection (parallel ray outside)" << endl;
  }

  cerr << "All intersect_2d tests passed!" << endl;
  return 0;
}
