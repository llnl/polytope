// 3D segment-face intersection on PLC unit test.

#include <iostream>
#include <vector>
#include <set>
#include <array>
#include <stdlib.h>
#include <limits>
#include <sstream>
#include <ctime>
#include <cmath>

#include "polytope.hh"
#include "intersect.hh"
#include "polytope_test_utilities.hh"

using namespace std;

//------------------------------------------------------------------------------
// Helper to check if a point is in the result set
//------------------------------------------------------------------------------
template<typename RealType>
bool containsPoint(const set<array<RealType, 3>>& result,
                   const RealType* point,
                   const RealType tol) {
  for (const auto& pt : result) {
    RealType dx = pt[0] - point[0];
    RealType dy = pt[1] - point[1];
    RealType dz = pt[2] - point[2];
    if (dx*dx + dy*dy + dz*dz < tol*tol) {
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

  // Create a cubic PLC with one cubic hole inside.
  // Outer cube: [-5, 5]^3
  // Inner cube (hole): [-1, 1]^3
  const unsigned numVertices = 16;
  double vertices[48] = {
    // Outer cube vertices (0-7)
    -5.0, -5.0, -5.0,  // 0
     5.0, -5.0, -5.0,  // 1
     5.0,  5.0, -5.0,  // 2
    -5.0,  5.0, -5.0,  // 3
    -5.0, -5.0,  5.0,  // 4
     5.0, -5.0,  5.0,  // 5
     5.0,  5.0,  5.0,  // 6
    -5.0,  5.0,  5.0,  // 7
    // Inner cube vertices (8-15) - hole
    -1.0, -1.0, -1.0,  // 8
     1.0, -1.0, -1.0,  // 9
     1.0,  1.0, -1.0,  // 10
    -1.0,  1.0, -1.0,  // 11
    -1.0, -1.0,  1.0,  // 12
     1.0, -1.0,  1.0,  // 13
     1.0,  1.0,  1.0,  // 14
    -1.0,  1.0,  1.0   // 15
  };

  polytope::PLC<3> plc;

  // Outer cube has 6 faces (each face is a quad)
  plc.facets.resize(6);
  plc.facets[0] = {0, 1, 2, 3}; // bottom (z = -5)
  plc.facets[1] = {4, 7, 6, 5}; // top (z = 5)
  plc.facets[2] = {0, 4, 5, 1}; // front (y = -5)
  plc.facets[3] = {2, 6, 7, 3}; // back (y = 5)
  plc.facets[4] = {0, 3, 7, 4}; // left (x = -5)
  plc.facets[5] = {1, 5, 6, 2}; // right (x = 5)

  // Inner cube (hole) has 6 faces
  plc.holes.resize(1);
  plc.holes[0].resize(6);
  plc.holes[0][0] = {8, 9, 10, 11};   // bottom (z = -1)
  plc.holes[0][1] = {12, 15, 14, 13}; // top (z = 1)
  plc.holes[0][2] = {8, 12, 13, 9};   // front (y = -1)
  plc.holes[0][3] = {10, 14, 15, 11}; // back (y = 1)
  plc.holes[0][4] = {8, 11, 15, 12};  // left (x = -1)
  plc.holes[0][5] = {9, 13, 14, 10};  // right (x = 1)

  { // Test 1: 1 intersection with outside boundary
    double p1[3] = {-10.0, 0.0, 0.0};
    double p2[3] = {-4.0, 0.0, 0.0};
    auto result = intersect(p1, p2, numVertices, vertices, plc);
    POLY_CHECK2(result.size() == 1, "Test 1: size=" << result.size());
    double expected[3] = {-5.0, 0.0, 0.0};
    POLY_CHECK2(containsPoint(result, expected, tol), "Test 1: missing expected intersection");
    cerr << "Test 1 passed: 1 intersection with outside boundary" << endl;
  }

  { // Test 2: 2 intersections - entry and exit through opposite faces
    double p1[3] = {-6.0, 2.0, 0.0};
    double p2[3] = {6.0, 2.0, 0.0};
    auto result = intersect(p1, p2, numVertices, vertices, plc);
    POLY_CHECK2(result.size() == 2, "Test 2: size=" << result.size());
    double expected1[3] = {-5.0, 2.0, 0.0};
    double expected2[3] = {5.0, 2.0, 0.0};
    POLY_CHECK2(containsPoint(result, expected1, tol), "Test 2: missing entry intersection");
    POLY_CHECK2(containsPoint(result, expected2, tol), "Test 2: missing exit intersection");
    cerr << "Test 2 passed: 2 intersections through opposite faces" << endl;
  }

  { // Test 3: Intersection with corner vertex
    double p1[3] = {4.0, 4.0, -6.0};
    double p2[3] = {6.0, 6.0, -4.0};
    auto result = intersect(p1, p2, numVertices, vertices, plc);
    POLY_CHECK2(result.size() == 1, "Test 3: size=" << result.size());
    double expected[3] = {5.0, 5.0, -5.0};
    POLY_CHECK2(containsPoint(result, expected, tol), "Test 3: missing corner intersection");
    cerr << "Test 3 passed: intersection at corner" << endl;
  }

  { // Test 4: Ray passing through hole (4 intersections total)
    double p1[3] = {0.0, 0.0, -6.0};
    double p2[3] = {0.0, 0.0, 6.0};
    auto result = intersect(p1, p2, numVertices, vertices, plc);
    POLY_CHECK2(result.size() == 4, "Test 4: size=" << result.size() << " (expected 4)");
    double expected[4][3] = {
      {0.0, 0.0, -5.0},  // outer bottom
      {0.0, 0.0, -1.0},  // hole bottom
      {0.0, 0.0, 1.0},   // hole top
      {0.0, 0.0, 5.0}    // outer top
    };
    for (int i = 0; i < 4; ++i) {
      POLY_CHECK2(containsPoint(result, expected[i], tol),
                  "Test 4: missing intersection " << i);
    }
    cerr << "Test 4 passed: 4 intersections through hole" << endl;
  }

  { // Test 5: Ray along edge
    double p1[3] = {5.0, 0.0, -6.0};
    double p2[3] = {5.0, 0.0, 6.0};
    auto result = intersect(p1, p2, numVertices, vertices, plc);
    POLY_CHECK2(result.size() == 2, "Test 5: size=" << result.size() << " (expected 2)");
    double expected1[3] = {5.0, 0.0, -5.0};
    double expected2[3] = {5.0, 0.0, 5.0};
    POLY_CHECK2(containsPoint(result, expected1, tol), "Test 5: missing bottom intersection");
    POLY_CHECK2(containsPoint(result, expected2, tol), "Test 5: missing top intersection");
    cerr << "Test 5 passed: 2 intersections along edge" << endl;
  }

  { // Test 6: Diagonal ray through cube (should hit 4 intersections with hole)
    double p1[3] = {-6.0, -6.0, -6.0};
    double p2[3] = {6.0, 6.0, 6.0};
    auto result = intersect(p1, p2, numVertices, vertices, plc);
    POLY_CHECK2(result.size() == 4, "Test 6: size=" << result.size() << " (expected 4)");
    double expected[4][3] = {
      {-5.0, -5.0, -5.0},  // outer entry
      {-1.0, -1.0, -1.0},  // hole entry
      {1.0, 1.0, 1.0},     // hole exit
      {5.0, 5.0, 5.0}      // outer exit
    };
    for (int i = 0; i < 4; ++i) {
      POLY_CHECK2(containsPoint(result, expected[i], tol),
                  "Test 6: missing diagonal intersection " << i);
    }
    cerr << "Test 6 passed: 4 diagonal intersections" << endl;
  }

  { // Test 7: Ray parallel to face (no intersection)
    double p1[3] = {-6.0, 0.0, 0.0};
    double p2[3] = {-6.0, 3.0, 3.0};
    auto result = intersect(p1, p2, numVertices, vertices, plc);
    POLY_CHECK2(result.size() == 0, "Test 7: size=" << result.size() << " (expected 0)");
    cerr << "Test 7 passed: no intersection (parallel ray outside)" << endl;
  }

  cerr << "All intersect_3d tests passed!" << endl;
  return 0;
}
