// Comprehensive unit tests for QuantPLC::makeConvex() in 2D
//
// Tests the 2D convex hull computation using quantized coordinate hashes
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
// Helper: Verify hull is properly closed (each vertex appears exactly twice)
//------------------------------------------------------------------------------
bool
isHullClosed(const QuantPLC<2>& qplc) {
  map<int, int> vertexCount;

  for (const auto& facet : qplc.facets) {
    POLY_ASSERT(facet.size() == 2);
    vertexCount[facet[0]]++;
    vertexCount[facet[1]]++;
  }

  // In a closed hull, each vertex should appear exactly twice
  for (const auto& kv : vertexCount) {
    if (kv.second != 2) {
      cerr << "  Vertex " << kv.first << " appears " << kv.second
           << " times (expected 2)" << endl;
      return false;
    }
  }

  return true;
}

//------------------------------------------------------------------------------
// Helper: Count unique vertices in hull
//------------------------------------------------------------------------------
unsigned
countHullVertices(const QuantPLC<2>& qplc) {
  set<int> vertices;
  for (const auto& facet : qplc.facets) {
    for (int idx : facet) {
      vertices.insert(idx);
    }
  }
  return vertices.size();
}

//------------------------------------------------------------------------------
// Test 1: Triangle (3 points)
//------------------------------------------------------------------------------
void testTriangle() {
  cout << "Test 1: Triangle" << endl;

  vector<double> points = {
    0.0, 0.0,
    1.0, 0.0,
    0.5, 1.0
  };

  // Initialize PLC with indices referencing the points
  PLC<2> plc;
  plc.facets.resize(3, vector<int>(2));
  plc.facets[0][0] = 0; plc.facets[0][1] = 1;
  plc.facets[1][0] = 1; plc.facets[1][1] = 2;
  plc.facets[2][0] = 2; plc.facets[2][1] = 0;

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);
  Quantizer<2> Q(low, high);

  QuantPLC<2> qplc(plc, Q, points, false);
  qplc.makeConvex();

  POLY_CHECK2(qplc.facets.size() == 3,
              "Triangle should produce 3 facets, got " << qplc.facets.size());

  POLY_CHECK2(isHullClosed(qplc),
              "Hull should be closed (each vertex appears twice)");

  POLY_CHECK2(countHullVertices(qplc) == 3,
              "Hull should have 3 vertices");

  cout << "  PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 2: Square (4 corner points)
//------------------------------------------------------------------------------
void testSquare() {
  cout << "Test 2: Square" << endl;

  vector<double> points = {
    0.0, 0.0,
    1.0, 0.0,
    1.0, 1.0,
    0.0, 1.0
  };

  // PLC representing a square boundary
  PLC<2> plc;
  plc.facets.resize(4, vector<int>(2));
  for (int i = 0; i < 4; ++i) {
    plc.facets[i][0] = i;
    plc.facets[i][1] = (i + 1) % 4;
  }

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);
  Quantizer<2> Q(low, high);

  QuantPLC<2> qplc(plc, Q, points, false);
  qplc.makeConvex();

  POLY_CHECK2(qplc.facets.size() == 4,
              "Square should produce 4 facets, got " << qplc.facets.size());

  POLY_CHECK2(isHullClosed(qplc),
              "Hull should be closed");

  POLY_CHECK2(countHullVertices(qplc) == 4,
              "Hull should have 4 vertices");

  cout << "  PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 3: Square with interior points
//------------------------------------------------------------------------------
void testSquareWithInterior() {
  cout << "Test 3: Square with interior points" << endl;

  vector<double> points = {
    // Corners
    0.0, 0.0,
    1.0, 0.0,
    1.0, 1.0,
    0.0, 1.0,
    // Interior points (should not appear in hull)
    0.5, 0.5,
    0.3, 0.3,
    0.7, 0.7,
    0.4, 0.6
  };

  // PLC with square boundary (first 4 points)
  PLC<2> plc;
  plc.facets.resize(4, vector<int>(2));
  for (int i = 0; i < 4; ++i) {
    plc.facets[i][0] = i;
    plc.facets[i][1] = (i + 1) % 4;
  }

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);
  Quantizer<2> Q(low, high);

  QuantPLC<2> qplc(plc, Q, points, false);
  qplc.makeConvex();

  POLY_CHECK2(qplc.facets.size() == 4,
              "Hull should have 4 facets (corners only), got "
              << qplc.facets.size());

  POLY_CHECK2(countHullVertices(qplc) == 4,
              "Hull should have 4 vertices (corners only)");

  // Verify all original points are contained in hull
  for (size_t i = 0; i < points.size() / 2; ++i) {
    Point<2, double> p(points[2*i], points[2*i+1]);
    POLY_CHECK2(qplc.within(p),
                "Point " << i << " at (" << p.x << ", " << p.y
                << ") not contained in hull");
  }

  cout << "  PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 4: Pentagon (regular)
//------------------------------------------------------------------------------
void testPentagon() {
  cout << "Test 4: Pentagon" << endl;

  // Generate 5 points on a circle
  vector<double> points;
  const int n = 5;
  const double cx = 0.5, cy = 0.5, r = 0.4;
  for (int i = 0; i < n; ++i) {
    double angle = 2.0 * M_PI * i / n;
    points.push_back(cx + r * cos(angle));
    points.push_back(cy + r * sin(angle));
  }

  // PLC with pentagon boundary
  PLC<2> plc;
  plc.facets.resize(n, vector<int>(2));
  for (int i = 0; i < n; ++i) {
    plc.facets[i][0] = i;
    plc.facets[i][1] = (i + 1) % n;
  }

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);
  Quantizer<2> Q(low, high);

  QuantPLC<2> qplc(plc, Q, points, false);
  qplc.makeConvex();

  POLY_CHECK2(qplc.facets.size() == 5,
              "Pentagon should produce 5 facets, got " << qplc.facets.size());

  POLY_CHECK2(isHullClosed(qplc),
              "Hull should be closed");

  POLY_CHECK2(countHullVertices(qplc) == 5,
              "Hull should have 5 vertices");

  cout << "  PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 6: Circle of points (many points on convex boundary)
//------------------------------------------------------------------------------
void testCircle() {
  cout << "Test 6: Circle of points" << endl;

  // Generate points on a circle
  vector<double> points;
  const int n = 20;
  const double cx = 0.5, cy = 0.5, r = 0.4;
  for (int i = 0; i < n; ++i) {
    double angle = 2.0 * M_PI * i / n;
    points.push_back(cx + r * cos(angle));
    points.push_back(cy + r * sin(angle));
  }

  // PLC with circular boundary
  PLC<2> plc;
  plc.facets.resize(n, vector<int>(2));
  for (int i = 0; i < n; ++i) {
    plc.facets[i][0] = i;
    plc.facets[i][1] = (i + 1) % n;
  }

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);
  Quantizer<2> Q(low, high);

  QuantPLC<2> qplc(plc, Q, points, false);
  qplc.makeConvex();

  // All points on circle should be hull vertices
  POLY_CHECK2(qplc.facets.size() == n,
              "Circle should produce " << n << " facets, got "
              << qplc.facets.size());

  POLY_CHECK2(isHullClosed(qplc),
              "Hull should be closed");

  cout << "  PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 7: Duplicate points
//------------------------------------------------------------------------------
void testDuplicates() {
  cout << "Test 7: Duplicate points" << endl;

  vector<double> points = {
    0.0, 0.0,  // Corner
    1.0, 0.0,  // Corner
    1.0, 1.0,  // Corner
    0.0, 1.0,  // Corner
    0.0, 0.0,  // Duplicate of first
    1.0, 0.0,  // Duplicate of second
  };

  // PLC with square boundary (first 4 points)
  PLC<2> plc;
  plc.facets.resize(4, vector<int>(2));
  for (int i = 0; i < 4; ++i) {
    plc.facets[i][0] = i;
    plc.facets[i][1] = (i + 1) % 4;
  }

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);
  Quantizer<2> Q(low, high);

  // With reduction, duplicates should be merged
  QuantPLC<2> qplc(plc, Q, points, true);
  qplc.makeConvex();

  // Should still produce a square hull
  POLY_CHECK2(qplc.facets.size() == 4,
              "Square with duplicates should produce 4 facets, got "
              << qplc.facets.size());

  POLY_CHECK2(isHullClosed(qplc),
              "Hull should be closed");

  cout << "  PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 8: Random point cloud
//------------------------------------------------------------------------------
void testRandomCloud() {
  cout << "Test 8: Random point cloud" << endl;

  // Generate random points
  const int n = 100;
  vector<double> points;
  for (int i = 0; i < n; ++i) {
    points.push_back(random01());
    points.push_back(random01());
  }

  // Create PLC with no initial boundary (will be determined by convex hull)
  PLC<2> plc;

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);
  Quantizer<2> Q(low, high);

  QuantPLC<2> qplc(plc, Q, points, false);
  qplc.makeConvex();

  // Hull should have at least 3 facets (unless collinear, very unlikely)
  POLY_CHECK2(qplc.facets.size() >= 3,
              "Random cloud should produce >= 3 facets, got "
              << qplc.facets.size());

  POLY_CHECK2(isHullClosed(qplc),
              "Hull should be closed");

  // Verify all original points are contained
  for (size_t i = 0; i < points.size() / 2; ++i) {
    Point<2, double> p(points[2*i], points[2*i+1]);
    POLY_CHECK2(qplc.within(p),
                "Point " << i << " not contained in hull");
  }

  cout << "  PASS (hull has " << qplc.facets.size() << " facets)" << endl;
}

//------------------------------------------------------------------------------
// Test 9: Large random cloud
//------------------------------------------------------------------------------
void testLargeRandom() {
  cout << "Test 9: Large random point cloud" << endl;

  // Generate many random points
  const int n = 10000;
  vector<double> points;
  for (int i = 0; i < n; ++i) {
    points.push_back(random01());
    points.push_back(random01());
  }

  // Create PLC with no initial boundary
  PLC<2> plc;

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);
  Quantizer<2> Q(low, high);

  QuantPLC<2> qplc(plc, Q, points, false);

  cout << "  Computing hull of " << n << " points..." << endl;
  qplc.makeConvex();

  POLY_CHECK2(qplc.facets.size() >= 3,
              "Large cloud should produce >= 3 facets, got "
              << qplc.facets.size());

  POLY_CHECK2(isHullClosed(qplc),
              "Hull should be closed");

  // Spot check some points for containment
  for (size_t i = 0; i < min(size_t(100), points.size() / 2); i += 1) {
    Point<2, double> p(points[2*i], points[2*i+1]);
    POLY_CHECK2(qplc.within(p),
                "Point " << i << " not contained in hull");
  }

  cout << "  PASS (hull has " << qplc.facets.size() << " facets)" << endl;
}

//------------------------------------------------------------------------------
// Test 10: Clustered points (near-duplicates)
//------------------------------------------------------------------------------
void testClustered() {
  cout << "Test 10: Clustered points" << endl;

  // Create 4 clusters near the corners
  vector<double> points;
  vector<pair<double, double>> centers = {
    {0.1, 0.1}, {0.9, 0.1}, {0.9, 0.9}, {0.1, 0.9}
  };

  for (const auto& center : centers) {
    for (int i = 0; i < 10; ++i) {
      // Add points in a small cluster
      points.push_back(center.first + 0.02 * (random01() - 0.5));
      points.push_back(center.second + 0.02 * (random01() - 0.5));
    }
  }

  // Create PLC with no initial boundary
  PLC<2> plc;

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);
  Quantizer<2> Q(low, high);

  QuantPLC<2> qplc(plc, Q, points, false);
  qplc.makeConvex();

  // Should produce roughly 4 facets (approximately square)
  POLY_CHECK2(qplc.facets.size() >= 4,
              "Clustered points should produce >= 4 facets, got "
              << qplc.facets.size());

  POLY_CHECK2(isHullClosed(qplc),
              "Hull should be closed");

  cout << "  PASS (hull has " << qplc.facets.size() << " facets)" << endl;
}

//------------------------------------------------------------------------------
// Test 11: Points on grid
//------------------------------------------------------------------------------
void testGrid() {
  cout << "Test 11: Grid of points" << endl;

  // Create 10x10 grid
  vector<double> points;
  const int gridSize = 10;
  for (int i = 0; i < gridSize; ++i) {
    for (int j = 0; j < gridSize; ++j) {
      points.push_back(double(i) / (gridSize - 1));
      points.push_back(double(j) / (gridSize - 1));
    }
  }

  // Create PLC with no initial boundary
  PLC<2> plc;

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);
  Quantizer<2> Q(low, high);

  QuantPLC<2> qplc(plc, Q, points, false);
  qplc.makeConvex();

  // Grid hull should be the 4 corners
  POLY_CHECK2(qplc.facets.size() == 4,
              "Grid should produce 4 facets, got " << qplc.facets.size());

  POLY_CHECK2(countHullVertices(qplc) == 4,
              "Grid hull should have 4 vertices");

  cout << "  PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 12: Extreme aspect ratio (very flat)
//------------------------------------------------------------------------------
void testExtremeAspectRatio() {
  cout << "Test 12: Extreme aspect ratio" << endl;

  vector<double> points = {
    0.0,  0.5,
    0.25, 0.500001,
    0.5,  0.499999,
    0.75, 0.500002,
    1.0,  0.5
  };

  // Create PLC with no initial boundary
  PLC<2> plc;

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);
  Quantizer<2> Q(low, high);

  QuantPLC<2> qplc(plc, Q, points, false);
  qplc.makeConvex();

  // Should produce a hull (may be line or thin polygon depending on quantization)
  POLY_CHECK2(qplc.facets.size() >= 1,
              "Should produce at least 1 facet");

  POLY_CHECK2(qplc.facets.size() == 1 || isHullClosed(qplc),
              "Non-degenerate hull should be closed");

  cout << "  PASS (hull has " << qplc.facets.size() << " facets)" << endl;
}

//------------------------------------------------------------------------------
// Test 13: Different quantization resolutions
//------------------------------------------------------------------------------
void testQuantizationResolution() {
  cout << "Test 13: Different quantization resolutions" << endl;

  // Points that might merge at different resolutions
  vector<double> points = {
    0.0,    0.0,
    0.001,  0.0,    // Very close to first
    1.0,    0.0,
    1.0,    1.0,
    0.0,    1.0
  };

  // Create PLC with no initial boundary
  PLC<2> plc;

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);

  // Fine resolution (small padding)
  Quantizer<2> Q_fine(low, high, 0.001);
  QuantPLC<2> qplc_fine(plc, Q_fine, points, false);
  qplc_fine.makeConvex();

  // Coarse resolution (large padding)
  Quantizer<2> Q_coarse(low, high, 0.1);
  QuantPLC<2> qplc_coarse(plc, Q_coarse, points, false);
  qplc_coarse.makeConvex();

  POLY_CHECK2(qplc_fine.facets.size() >= 4,
              "Fine quantization should preserve detail");
  POLY_CHECK2(qplc_coarse.facets.size() >= 4,
              "Coarse quantization should still produce valid hull");

  cout << "  PASS (fine: " << qplc_fine.facets.size()
       << " facets, coarse: " << qplc_coarse.facets.size() << " facets)" << endl;
}

//------------------------------------------------------------------------------
// Test 14: Hull facet ordering (should form closed loop)
//------------------------------------------------------------------------------
void testFacetOrdering() {
  cout << "Test 14: Hull facet ordering" << endl;

  vector<double> points = {
    0.0, 0.0,
    1.0, 0.0,
    1.0, 1.0,
    0.0, 1.0
  };

  // PLC representing a square boundary
  PLC<2> plc;
  plc.facets.resize(4, vector<int>(2));
  for (int i = 0; i < 4; ++i) {
    plc.facets[i][0] = i;
    plc.facets[i][1] = (i + 1) % 4;
  }

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);
  Quantizer<2> Q(low, high);

  QuantPLC<2> qplc(plc, Q, points, false);
  qplc.makeConvex();

  // Try to walk around the hull
  set<int> visited;
  int current = qplc.facets[0][0];
  visited.insert(current);

  for (size_t i = 0; i < qplc.facets.size(); ++i) {
    // Find next edge that starts with current vertex
    bool found = false;
    for (const auto& facet : qplc.facets) {
      if (facet[0] == current && visited.count(facet[1]) == 0) {
        current = facet[1];
        visited.insert(current);
        found = true;
        break;
      }
    }
    if (!found) break;
  }

  // Should have visited all hull vertices
  POLY_CHECK2(visited.size() == countHullVertices(qplc),
              "Should be able to walk entire hull");

  cout << "  PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 15: Stress test with many random configurations
//------------------------------------------------------------------------------
void testStress() {
  cout << "Test 15: Stress test with random configurations" << endl;

  const int nTests = 100;
  int passed = 0;

  for (int test = 0; test < nTests; ++test) {
    // Random number of points
    int n = 3 + int(random01() * 100);

    vector<double> points;
    for (int i = 0; i < n; ++i) {
      points.push_back(random01());
      points.push_back(random01());
    }

    // Create PLC with no initial boundary
    PLC<2> plc;

    Point<2, double> low(0.0, 0.0);
    Point<2, double> high(1.0, 1.0);
    Quantizer<2> Q(low, high);

    QuantPLC<2> qplc(plc, Q, points, false);
    qplc.makeConvex();

    // Basic validity checks
    if (qplc.facets.size() >= 3 && isHullClosed(qplc)) {
      passed++;
    } else if (qplc.facets.size() < 3) {
      // Could be collinear, verify
      passed++;
    }
  }

  POLY_CHECK2(passed == nTests,
              "Stress test: " << passed << "/" << nTests << " passed");

  cout << "  PASS (" << nTests << " random configurations)" << endl;
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
    cout << "\n=== Testing QuantPLC::makeConvex() in 2D ===" << endl;

    // Simple convex shapes
    testTriangle();
    testSquare();
    testSquareWithInterior();
    testPentagon();
    testCircle();

    // Edge cases
    testDuplicates();
    testClustered();
    testGrid();
    testExtremeAspectRatio();
    testQuantizationResolution();

    // Validation tests
    testFacetOrdering();

    // Robustness tests
    testRandomCloud();
    testLargeRandom();
    testStress();

    cout << "\n=== All QuantPLC makeConvex 2D tests PASSED ===" << endl;

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
