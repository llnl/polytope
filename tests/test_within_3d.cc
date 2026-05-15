// 3D within on PLC unit test.

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <limits>
#include <sstream>
#include <ctime>

#include "polytope.hh"
#include "within.hh"
#include "polytope_test_utilities.hh"

using namespace std;

//------------------------------------------------------------------------------
// The test itself.
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

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

  unsigned testNum = 1;

  { // Test 1: Inside - fully interior
    double p[3] = {-4.0, 1.0, 2.0};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << testNum << ": inside=" << inside << " expected=" << answer);
    cerr << "Test " << testNum << " passed: fully interior point" << endl;
    ++testNum;
  }

  { // Test 2: Outside - completely outside outer boundary
    double p[3] = {-6.0, -3.0, 2.0};
    bool answer = false;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << testNum);
    cerr << "Test " << testNum << " passed: fully outside point" << endl;
    ++testNum;
  }

  { // Test 3: Outside - within hole
    double p[3] = {0.0, 0.0, 0.0};
    bool answer = false;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << testNum);
    cerr << "Test " << testNum << " passed: point within hole" << endl;
    ++testNum;
  }

  { // Test 4: Inside - on outer boundary (face)
    double p[3] = {-5.0, 0.0, 0.0};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << testNum);
    cerr << "Test " << testNum << " passed: on outer boundary face" << endl;
    ++testNum;
  }

  { // Test 5: Inside - on outer boundary (edge)
    double p[3] = {-5.0, -5.0, 0.0};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << testNum);
    cerr << "Test " << testNum << " passed: on outer boundary edge" << endl;
    ++testNum;
  }

  { // Test 6: Inside - on outer boundary (vertex)
    double p[3] = {-5.0, -5.0, -5.0};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << testNum);
    cerr << "Test " << testNum << " passed: on outer boundary vertex" << endl;
    ++testNum;
  }

  { // Test 7: Inside - on hole boundary (convention: on hole boundary = inside)
    double p[3] = {-1.0, 0.0, 0.0};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << testNum);
    cerr << "Test " << testNum << " passed: on hole boundary" << endl;
    ++testNum;
  }

  { // Test 8: Inside - between outer boundary and hole
    double p[3] = {-3.0, 0.0, 0.0};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << testNum);
    cerr << "Test " << testNum << " passed: between boundary and hole" << endl;
    ++testNum;
  }

  { // Test 9: Inside - near corner of outer boundary
    double p[3] = {4.5, 4.5, 4.5};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << testNum);
    cerr << "Test " << testNum << " passed: near outer corner" << endl;
    ++testNum;
  }

  { // Test 10: Outside - beyond each axis
    double p[3] = {6.0, 0.0, 0.0};
    bool answer = false;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << testNum << " (x-axis)");
    cerr << "Test " << testNum << " passed: outside along x-axis" << endl;
    ++testNum;
  }

  { // Test 11: Outside - beyond y-axis
    double p[3] = {0.0, 6.0, 0.0};
    bool answer = false;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << testNum << " (y-axis)");
    cerr << "Test " << testNum << " passed: outside along y-axis" << endl;
    ++testNum;
  }

  { // Test 12: Outside - beyond z-axis
    double p[3] = {0.0, 0.0, 6.0};
    bool answer = false;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << testNum << " (z-axis)");
    cerr << "Test " << testNum << " passed: outside along z-axis" << endl;
    ++testNum;
  }

  { // Test 13: Test a tetrahedron (non-cube geometry)
    // Create a simple tetrahedron
    const unsigned nTetVerts = 4;
    double tetVerts[12] = {
      0.0, 0.0, 0.0,   // 0
      1.0, 0.0, 0.0,   // 1
      0.5, 1.0, 0.0,   // 2
      0.5, 0.5, 1.0    // 3
    };

    polytope::PLC<3> tetPLC;
    tetPLC.facets.resize(4);
    tetPLC.facets[0] = {0, 2, 1};  // base
    tetPLC.facets[1] = {0, 1, 3};  // side 1
    tetPLC.facets[2] = {0, 3, 2};  // side 2
    tetPLC.facets[3] = {1, 2, 3};  // side 3

    { // Inside tetrahedron
      double p[3] = {0.4, 0.3, 0.2};
      bool answer = true;
      const bool inside = within(p, nTetVerts, tetVerts, tetPLC);
      POLY_CHECK2(inside == answer, "Test " << testNum << " (inside tet)");
    }

    { // Outside tetrahedron
      double p[3] = {2.0, 0.0, 0.0};
      bool answer = false;
      const bool inside = within(p, nTetVerts, tetVerts, tetPLC);
      POLY_CHECK2(inside == answer, "Test " << testNum << " (outside tet)");
    }

    { // On tetrahedron face
      double p[3] = {0.5, 0.5, 0.0};
      bool answer = true;
      const bool inside = within(p, nTetVerts, tetVerts, tetPLC);
      POLY_CHECK2(inside == answer, "Test " << testNum << " (on tet face)");
    }

    cerr << "Test " << testNum << " passed: tetrahedron tests" << endl;
    ++testNum;
  }

  cerr << "All within_3d tests passed!" << endl;
  return 0;
}
