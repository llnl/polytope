// Comprehensive unit tests for QuantPLC::makeConvex() in 3D
//
// Tests the 3D convex hull computation using quantized coordinate hashes
// with the QuantPLC class.

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <cassert>
#include <cmath>
#include <algorithm>

#include "polytope.hh"
#include "QuantPLC.hh"
#include "PLC.hh"
#include "Point.hh"
#include "HashKey.hh"
#include "Quantizer.hh"
#include "polytope_test_utilities.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

namespace {

//------------------------------------------------------------------------------
// Helper: Count unique vertices in hull
//------------------------------------------------------------------------------
unsigned
countHullVertices3D(const QuantPLC<3>& qplc) {
  set<int> vertices;
  for (const auto& facet : qplc.facets) {
    for (int idx : facet) {
      vertices.insert(idx);
    }
  }
  return vertices.size();
}

//------------------------------------------------------------------------------
// Helper: Check hull is properly formed (Euler characteristic)
//------------------------------------------------------------------------------
bool
checkEulerCharacteristic(const QuantPLC<3>& qplc) {
  // For a closed convex polyhedron: V - E + F = 2
  unsigned V = countHullVertices3D(qplc);
  unsigned F = qplc.facets.size();

  // Count edges (each edge shared by exactly 2 facets)
  map<pair<int, int>, int> edgeCount;
  for (const auto& facet : qplc.facets) {
    POLY_ASSERT(facet.size() == 3);
    for (int i = 0; i < 3; ++i) {
      int a = facet[i];
      int b = facet[(i+1)%3];
      auto edge = make_pair(min(a, b), max(a, b));
      edgeCount[edge]++;
    }
  }
  unsigned E = edgeCount.size();

  // Check each edge is shared by exactly 2 facets
  for (const auto& kv : edgeCount) {
    if (kv.second != 2) {
      cerr << "  Edge (" << kv.first.first << "," << kv.first.second
           << ") appears " << kv.second << " times (expected 2)" << endl;
      return false;
    }
  }

  // Check Euler characteristic
  if (V - E + F != 2) {
    cerr << "  Euler characteristic: V=" << V << " E=" << E << " F=" << F
         << " V-E+F=" << (V-E+F) << " (expected 2)" << endl;
    return false;
  }

  return true;
}

//------------------------------------------------------------------------------
// Test 1: Tetrahedron (4 points)
//------------------------------------------------------------------------------
void testTetrahedron() {
  cout << "Test 1: Tetrahedron" << endl;

  vector<double> points = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.5, 1.0, 0.0,
    0.5, 0.5, 1.0
  };

  PLC<3> plc;
  plc.facets.resize(4, vector<int>(3));
  plc.facets[0] = {0, 1, 2};
  plc.facets[1] = {0, 1, 3};
  plc.facets[2] = {1, 2, 3};
  plc.facets[3] = {2, 0, 3};

  Point<3, double> low(0.0, 0.0, 0.0);
  Point<3, double> high(1.0, 1.0, 1.0);
  Quantizer<3> Q(low, high);

  QuantPLC<3> qplc(plc, Q, points, false);
  qplc.makeConvex();

  POLY_CHECK2(qplc.facets.size() == 4,
              "Tetrahedron should have 4 facets, got " << qplc.facets.size());

  POLY_CHECK2(countHullVertices3D(qplc) == 4,
              "Tetrahedron should have 4 vertices");

  POLY_CHECK2(checkEulerCharacteristic(qplc),
              "Hull should satisfy Euler characteristic");

  cout << "  PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 5: Cube (8 corner points)
//------------------------------------------------------------------------------
void testCube() {
  cout << "Test 5: Cube" << endl;

  vector<double> points = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    1.0, 0.0, 1.0,
    1.0, 1.0, 1.0,
    0.0, 1.0, 1.0
  };

  // Create PLC with cube faces
  PLC<3> plc;
  plc.facets = {
    {0, 1, 2}, {0, 2, 3},  // Bottom
    {4, 5, 6}, {4, 6, 7},  // Top
    {0, 1, 5}, {0, 5, 4},  // Front
    {2, 3, 7}, {2, 7, 6},  // Back
    {0, 3, 7}, {0, 7, 4},  // Left
    {1, 2, 6}, {1, 6, 5}   // Right
  };

  Point<3, double> low(0.0, 0.0, 0.0);
  Point<3, double> high(1.0, 1.0, 1.0);
  Quantizer<3> Q(low, high);

  QuantPLC<3> qplc(plc, Q, points, false);
  qplc.makeConvex();

  POLY_CHECK2(qplc.facets.size() == 12,
              "Cube should have 12 triangular facets, got " << qplc.facets.size());

  POLY_CHECK2(countHullVertices3D(qplc) == 8,
              "Cube should have 8 vertices");

  POLY_CHECK2(checkEulerCharacteristic(qplc),
              "Hull should satisfy Euler characteristic");

  cout << "  PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 3: Cube with interior points
//------------------------------------------------------------------------------
void testCubeWithInterior() {
  cout << "Test 3: Cube with interior points" << endl;

  vector<double> points = {
    // Corners
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    1.0, 0.0, 1.0,
    1.0, 1.0, 1.0,
    0.0, 1.0, 1.0,
    // Interior points
    0.5, 0.5, 0.5,
    0.3, 0.3, 0.3,
    0.7, 0.7, 0.7
  };

  PLC<3> plc;
  plc.facets = {
    {0, 1, 2}, {0, 2, 3},
    {4, 5, 6}, {4, 6, 7},
    {0, 1, 5}, {0, 5, 4},
    {2, 3, 7}, {2, 7, 6},
    {0, 3, 7}, {0, 7, 4},
    {1, 2, 6}, {1, 6, 5}
  };

  Point<3, double> low(0.0, 0.0, 0.0);
  Point<3, double> high(1.0, 1.0, 1.0);
  Quantizer<3> Q(low, high);

  QuantPLC<3> qplc(plc, Q, points, false);
  qplc.makeConvex();

  POLY_CHECK2(qplc.facets.size() == 12,
              "Cube hull should have 12 facets (corners only)");

  POLY_CHECK2(countHullVertices3D(qplc) == 8,
              "Cube hull should have 8 vertices (corners only)");

  // Verify all original points are contained
  for (size_t i = 0; i < points.size() / 3; ++i) {
    Point<3, double> p(points[3*i], points[3*i+1], points[3*i+2]);
    POLY_CHECK2(qplc.within(p), "Point " << i << " not contained in hull");
  }

  cout << "  PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 4: Octahedron
//------------------------------------------------------------------------------
void testOctahedron() {
  cout << "Test 4: Octahedron" << endl;

  // 6 vertices at (+/-1, 0, 0), (0, +/-1, 0), (0, 0, +/-1)
  vector<double> points = {
    1.0, 0.0, 0.0,
    -1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, -1.0, 0.0,
    0.0, 0.0, 1.0,
    0.0, 0.0, -1.0
  };

  PLC<3> plc;
  // Octahedron has 8 triangular faces
  plc.facets = {
    {0, 2, 4}, {0, 4, 3}, {0, 3, 5}, {0, 5, 2},
    {1, 2, 4}, {1, 4, 3}, {1, 3, 5}, {1, 5, 2}
  };

  Point<3, double> low(-1.0, -1.0, -1.0);
  Point<3, double> high(1.0, 1.0, 1.0);
  Quantizer<3> Q(low, high);

  QuantPLC<3> qplc(plc, Q, points, false);
  qplc.makeConvex();

  POLY_CHECK2(qplc.facets.size() == 8,
              "Octahedron should have 8 facets");

  POLY_CHECK2(countHullVertices3D(qplc) == 6,
              "Octahedron should have 6 vertices");

  POLY_CHECK2(checkEulerCharacteristic(qplc),
              "Hull should satisfy Euler characteristic");

  cout << "  PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 5: Sphere approximation (icosahedron-like)
//------------------------------------------------------------------------------
void testSphere() {
  cout << "Test 5: Sphere approximation" << endl;

  // Generate points on a sphere
  vector<double> points;
  const int nTheta = 6;
  const int nPhi = 6;
  const double cx = 0.5, cy = 0.5, cz = 0.5, r = 0.4;

  for (int i = 0; i < nTheta; ++i) {
    for (int j = 0; j < nPhi; ++j) {
      double theta = M_PI * i / (nTheta - 1);
      double phi = 2.0 * M_PI * j / nPhi;
      points.push_back(cx + r * sin(theta) * cos(phi));
      points.push_back(cy + r * sin(theta) * sin(phi));
      points.push_back(cz + r * cos(theta));
    }
  }

  PLC<3> plc;

  Point<3, double> low(0.0, 0.0, 0.0);
  Point<3, double> high(1.0, 1.0, 1.0);
  Quantizer<3> Q(low, high);

  QuantPLC<3> qplc(plc, Q, points, false);
  qplc.makeConvex();

  POLY_CHECK2(qplc.facets.size() > 0,
              "Sphere approximation should have facets");

  POLY_CHECK2(checkEulerCharacteristic(qplc),
              "Hull should satisfy Euler characteristic");

  cout << "  PASS (hull has " << qplc.facets.size() << " facets)" << endl;
}

//------------------------------------------------------------------------------
// Test 6: Random point cloud
//------------------------------------------------------------------------------
void testRandomCloud() {
  cout << "Test 6: Random point cloud" << endl;

  // Generate random points
  const int n = 50;
  vector<double> points;
  for (int i = 0; i < n; ++i) {
    points.push_back(random01());
    points.push_back(random01());
    points.push_back(random01());
  }

  PLC<3> plc;

  Point<3, double> low(0.0, 0.0, 0.0);
  Point<3, double> high(1.0, 1.0, 1.0);
  Quantizer<3> Q(low, high);

  QuantPLC<3> qplc(plc, Q, points, false);
  qplc.makeConvex();

  POLY_CHECK2(qplc.facets.size() >= 4,
              "Random cloud should produce >= 4 facets");

  POLY_CHECK2(checkEulerCharacteristic(qplc),
              "Hull should satisfy Euler characteristic");

  // Verify all points are contained
  for (size_t i = 0; i < points.size() / 3; ++i) {
    Point<3, double> p(points[3*i], points[3*i+1], points[3*i+2]);
    POLY_CHECK2(qplc.within(p), "Point " << i << " not contained in hull");
  }

  cout << "  PASS (hull has " << qplc.facets.size() << " facets, "
       << countHullVertices3D(qplc) << " vertices)" << endl;
}

//------------------------------------------------------------------------------
// Test 7: Large random cloud
//------------------------------------------------------------------------------
void testLargeRandom() {
  cout << "Test 7: Large random point cloud" << endl;

  // Generate many random points
  const int n = 1000;
  vector<double> points;
  for (int i = 0; i < n; ++i) {
    points.push_back(random01());
    points.push_back(random01());
    points.push_back(random01());
  }

  PLC<3> plc;

  Point<3, double> low(0.0, 0.0, 0.0);
  Point<3, double> high(1.0, 1.0, 1.0);
  Quantizer<3> Q(low, high);

  QuantPLC<3> qplc(plc, Q, points, false);

  cout << "  Computing hull of " << n << " points..." << endl;
  qplc.makeConvex();

  POLY_CHECK2(qplc.facets.size() >= 4,
              "Large cloud should produce >= 4 facets");

  POLY_CHECK2(checkEulerCharacteristic(qplc),
              "Hull should satisfy Euler characteristic");

  // Spot check some points
  for (size_t i = 0; i < min(size_t(50), points.size() / 3); i += 1) {
    Point<3, double> p(points[3*i], points[3*i+1], points[3*i+2]);
    POLY_CHECK2(qplc.within(p), "Point " << i << " not contained in hull");
  }

  cout << "  PASS (hull has " << qplc.facets.size() << " facets, "
       << countHullVertices3D(qplc) << " vertices)" << endl;
}

//------------------------------------------------------------------------------
// Test 8: Clustered points
//------------------------------------------------------------------------------
void testClustered() {
  cout << "Test 8: Clustered points" << endl;

  // Create 8 clusters near cube corners
  vector<double> points;
  vector<tuple<double, double, double>> centers = {
    {0.1, 0.1, 0.1}, {0.9, 0.1, 0.1},
    {0.9, 0.9, 0.1}, {0.1, 0.9, 0.1},
    {0.1, 0.1, 0.9}, {0.9, 0.1, 0.9},
    {0.9, 0.9, 0.9}, {0.1, 0.9, 0.9}
  };

  for (const auto& center : centers) {
    for (int i = 0; i < 5; ++i) {
      points.push_back(get<0>(center) + 0.01 * (random01() - 0.5));
      points.push_back(get<1>(center) + 0.01 * (random01() - 0.5));
      points.push_back(get<2>(center) + 0.01 * (random01() - 0.5));
    }
  }

  PLC<3> plc;

  Point<3, double> low(0.0, 0.0, 0.0);
  Point<3, double> high(1.0, 1.0, 1.0);
  Quantizer<3> Q(low, high);

  QuantPLC<3> qplc(plc, Q, points, false);
  qplc.makeConvex();

  // Should produce approximately a cube hull
  POLY_CHECK2(qplc.facets.size() >= 8,
              "Clustered points should produce >= 8 facets");

  POLY_CHECK2(checkEulerCharacteristic(qplc),
              "Hull should satisfy Euler characteristic");

  cout << "  PASS (hull has " << qplc.facets.size() << " facets)" << endl;
}

//------------------------------------------------------------------------------
// Test 9: Grid of points
//------------------------------------------------------------------------------
void testGrid() {
  cout << "Test 9: Grid of points" << endl;

  // Create 5x5x5 grid
  vector<double> points;
  const int gridSize = 5;
  for (int i = 0; i < gridSize; ++i) {
    for (int j = 0; j < gridSize; ++j) {
      for (int k = 0; k < gridSize; ++k) {
        points.push_back(double(i) / (gridSize - 1));
        points.push_back(double(j) / (gridSize - 1));
        points.push_back(double(k) / (gridSize - 1));
      }
    }
  }

  PLC<3> plc;

  Point<3, double> low(0.0, 0.0, 0.0);
  Point<3, double> high(1.0, 1.0, 1.0);
  Quantizer<3> Q(low, high);

  QuantPLC<3> qplc(plc, Q, points, false);
  qplc.makeConvex();

  // Grid hull should be the 8 corners (cube)
  POLY_CHECK2(countHullVertices3D(qplc) == 8,
              "Grid hull should have 8 vertices");

  POLY_CHECK2(qplc.facets.size() == 12,
              "Grid hull (cube) should have 12 facets");

  cout << "  PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 10: Stress test
//------------------------------------------------------------------------------
void testStress() {
  cout << "Test 10: Stress test with random configurations" << endl;

  const int nTests = 50;
  int passed = 0;

  for (int test = 0; test < nTests; ++test) {
    int n = 10 + int(random01() * 50);

    vector<double> points;
    for (int i = 0; i < n; ++i) {
      points.push_back(random01());
      points.push_back(random01());
      points.push_back(random01());
    }

    PLC<3> plc;

    Point<3, double> low(0.0, 0.0, 0.0);
    Point<3, double> high(1.0, 1.0, 1.0);
    Quantizer<3> Q(low, high);

    QuantPLC<3> qplc(plc, Q, points, false);
    qplc.makeConvex();

    if (qplc.facets.size() >= 4 && checkEulerCharacteristic(qplc)) {
      passed++;
    }
  }

  POLY_CHECK2(passed == nTests,
              "Stress test: " << passed << "/" << nTests << " passed");

  cout << "  PASS (" << nTests << " random configurations)" << endl;
}

//------------------------------------------------------------------------------
// Test 11: QHull reproducibility with quantized coordinates
//------------------------------------------------------------------------------
void testReproducibility() {
  cout << "Test 11: QHull reproducibility with quantized coordinates" << endl;

  // Create a set of points for testing
  vector<double> points = {
    0.1, 0.2, 0.3,
    0.9, 0.1, 0.2,
    0.8, 0.9, 0.1,
    0.2, 0.8, 0.2,
    0.3, 0.3, 0.9,
    0.7, 0.2, 0.8,
    0.6, 0.7, 0.7,
    0.2, 0.6, 0.8,
    0.5, 0.5, 0.5,  // Interior point
    0.4, 0.4, 0.6,
    0.6, 0.3, 0.4
  };

  Point<3, double> low(0.0, 0.0, 0.0);
  Point<3, double> high(1.0, 1.0, 1.0);
  Quantizer<3> Q(low, high);

  // Run makeConvex multiple times and capture results
  const int nTrials = 10;
  vector<QuantPLC<3>> allQPLC;

  for (int trial = 0; trial < nTrials; ++trial) {
    PLC<3> plc;
    QuantPLC<3> qplc(plc, Q, points, false);
    qplc.makeConvex();

    allQPLC.push_back(qplc);
  }

  // Verify all trials produced identical results
  for (int trial = 1; trial < nTrials; ++trial) {
    POLY_CHECK2(QuantPLC<3>::compareFacets(allQPLC[trial], allQPLC[0]),
                "Trial " << trial << " produced different facets than trial 0");
  }

  // Verify Euler characteristic on first result
  PLC<3> plc;
  QuantPLC<3> qplc(plc, Q, points, false);
  qplc.makeConvex();
  POLY_CHECK2(checkEulerCharacteristic(qplc),
              "Hull should satisfy Euler characteristic");

  cout << "  PASS (" << nTrials << " trials: facets, and vertices all identical)" << endl;
}

//------------------------------------------------------------------------------
// Test 12: Translation invariance
//------------------------------------------------------------------------------
void testTranslation() {
  cout << "Test 12: Translation invariance" << endl;

  // Original points
  vector<double> points1 = {
    0.1, 0.2, 0.3,
    0.9, 0.1, 0.2,
    0.8, 0.9, 0.1,
    0.2, 0.8, 0.2,
    0.3, 0.3, 0.9,
    0.7, 0.2, 0.8,
    0.6, 0.7, 0.7,
    0.2, 0.6, 0.8,    
    0.5, 0.5, 0.5
  };

  // Amount to shift it by
  const double dx = 0.00001;
  double lo_point = 0.0;
  double hi_point = 0.95;
  const int nTrials = 10000;
  std::vector<QuantPLC<3>> allQPLC;

  for (int trial = 0; trial < nTrials; ++trial) {
    Point<3, double> low(lo_point, lo_point, lo_point);
    Point<3, double> hi(hi_point, hi_point, hi_point);
    PLC<3> plc;
    Quantizer<3> Q(low, hi);
    QuantPLC<3> qplc(plc, Q, points1);
    qplc.makeConvex();
    allQPLC.push_back(qplc);
    for (auto& v : points1) {
      v += dx;
    }
    lo_point += dx;
    hi_point += dx;
  }

  for (int trial = 1; trial < nTrials; ++trial) {
    POLY_CHECK2(QuantPLC<3>::compareFacets(allQPLC[trial], allQPLC[0]),
                "Trial " << trial << " produced a different facet set than trial 0");
  }

  cout << "  PASS (translated hull matches original hull)" << endl;
}

//------------------------------------------------------------------------------
// Test 13: Shuffle invariance
//------------------------------------------------------------------------------
void testShuffle() {
  cout << "Test 13: Shuffle invariance" << endl;

  // Original ordered points
  vector<double> points1 = {
    0.1, 0.2, 0.3,   // 0
    0.9, 0.1, 0.2,   // 1
    0.8, 0.9, 0.1,   // 2
    0.2, 0.8, 0.2,   // 3
    0.3, 0.3, 0.9,   // 4
    0.7, 0.2, 0.8,   // 5
    0.6, 0.7, 0.7,   // 6
    0.2, 0.6, 0.8,   // 7
    0.5, 0.5, 0.5    // 8
  };

  // Shuffled points with mapping
  vector<int> shuffle = {4, 1, 7, 0, 8, 3, 6, 2, 5};  // New order
  vector<int> invShuffle(9);  // invShuffle[new_idx] = old_idx
  for (size_t i = 0; i < shuffle.size(); ++i) {
    invShuffle[i] = shuffle[i];
  }

  vector<double> points2;
  for (int idx : shuffle) {
    points2.push_back(points1[3*idx + 0]);
    points2.push_back(points1[3*idx + 1]);
    points2.push_back(points1[3*idx + 2]);
  }

  Point<3, double> low(0.0, 0.0, 0.0);
  Point<3, double> high(1.0, 1.0, 1.0);
  Quantizer<3> Q(low, high);

  // Compute hulls
  PLC<3> plc1, plc2;
  QuantPLC<3> qplc1(plc1, Q, points1, false);
  QuantPLC<3> qplc2(plc2, Q, points2, false);
  qplc1.makeConvex();
  qplc2.makeConvex();

  POLY_CHECK2(QuantPLC<3>::compareFacets(qplc1, qplc2),
              "Shuffled hull should have same coordinate hashes as original");

  // cout << "  PASS (both hulls have " << verts1.size()
  //      << " vertices, " << qplc1.facets.size() << " facets)" << endl;
}

//------------------------------------------------------------------------------
// Test 14: Random reordering invariance
//------------------------------------------------------------------------------
void testReordering() {
  cout << "Test 14: Random reordering invariance" << endl;

  // Original points
  vector<double> points1 = {
    0.15, 0.25, 0.35,
    0.85, 0.15, 0.25,
    0.75, 0.85, 0.15,
    0.25, 0.75, 0.25,
    0.35, 0.35, 0.85,
    0.65, 0.25, 0.75,
    0.55, 0.65, 0.65,
    0.25, 0.55, 0.75,
    0.45, 0.45, 0.45,
    0.35, 0.65, 0.55
  };

  const int nPoints = points1.size() / 3;

  Point<3, double> low(0.0, 0.0, 0.0);
  Point<3, double> high(1.0, 1.0, 1.0);
  Quantizer<3> Q(low, high);

  // Compute reference hull
  PLC<3> plc1;
  QuantPLC<3> qplc1(plc1, Q, points1, false);
  qplc1.makeConvex();

  set<int> verts1;
  for (const auto& facet : qplc1.facets) {
    for (int idx : facet) verts1.insert(idx);
  }

  // Test multiple random reorderings
  const int nTrials = 5;
  for (int trial = 0; trial < nTrials; ++trial) {
    // Create random permutation
    vector<int> perm(nPoints);
    for (int i = 0; i < nPoints; ++i) perm[i] = i;
    for (int i = nPoints - 1; i > 0; --i) {
      int j = int(random01() * (i + 1));
      swap(perm[i], perm[j]);
    }

    // Inverse permutation for mapping back
    vector<int> invPerm(nPoints);
    for (int i = 0; i < nPoints; ++i) {
      invPerm[perm[i]] = i;
    }

    // Reorder points
    vector<double> points2;
    for (int idx : perm) {
      points2.push_back(points1[3*idx + 0]);
      points2.push_back(points1[3*idx + 1]);
      points2.push_back(points1[3*idx + 2]);
    }

    // Compute hull
    PLC<3> plc2;
    QuantPLC<3> qplc2(plc2, Q, points2, false);
    qplc2.makeConvex();

    POLY_CHECK2(QuantPLC<3>::compareFacets(qplc1, qplc2),
                "Reordering trial " << trial << " produced different hull vertices");

    POLY_CHECK2(qplc1.facets.size() == qplc2.facets.size(),
                "Reordering trial " << trial << " produced different facet count");
  }

  cout << "  PASS (" << nTrials << " random reorderings: "
       << verts1.size() << " vertices, "
       << qplc1.facets.size() << " facets, all consistent)" << endl;
}

} // anonymous namespace

//------------------------------------------------------------------------------
// Main
//------------------------------------------------------------------------------
int main(int argc, char** argv) {
#ifdef POLYTOPE_ENABLE_MPI
  MPI_Init(&argc, &argv);
#else
  POLY_CONTRACT_VAR(argc);
  POLY_CONTRACT_VAR(argv);
#endif

  try {
    cout << "\n=== Testing QuantPLC::makeConvex() in 3D ===" << endl;

    // Simple convex shapes
    testTetrahedron();
    testCube();
    testCubeWithInterior();
    testOctahedron();
    testSphere();

    // Point clouds
    testRandomCloud();
    testLargeRandom();
    testClustered();
    testGrid();

    // Robustness
    testStress();

    // Consistency and invariance
    testReproducibility();
    testTranslation();
    testShuffle();
    testReordering();

    cout << "\n=== All QuantPLC makeConvex 3D tests PASSED ===" << endl;

  } catch (const char* str) {
    cerr << "FAILED: " << str << endl;
#ifdef POLYTOPE_ENABLE_MPI
    MPI_Finalize();
#endif
    return 1;
  } catch (const std::exception& e) {
    cerr << "FAILED: " << e.what() << endl;
#ifdef POLYTOPE_ENABLE_MPI
    MPI_Finalize();
#endif
    return 1;
  }

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
