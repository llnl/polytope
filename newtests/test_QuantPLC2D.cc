// Comprehensive unit tests for QuantPLC<2>
//
// Tests the quantized PLC functionality for 2D including:
//   - Construction from real-space coordinates
//   - Quantization and deduplication
//   - Convex hull computation (2D)
//   - Facet (edge) ordering
//   - Collinear vertex removal
//   - Point containment (within) tests
//   - Intersection tests
//   - Reduction (removing unused points)
//   - Hash-based point comparison

#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <cassert>

#include "polytope.hh"
#include "QuantPLC.hh"
#include "Quantizer.hh"
#include "Intersections.hh"
#include "GeomUtils.hh"
#include "Point.hh"
#include "polytope_test_utilities.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace polytope;

namespace {

using RealType = double;
using RealPoint = Point2<RealType>;
using PLC = PLC<2>;
using QuantPLC2D = QuantPLC<2>;
using Quantizer2D = Quantizer<2>;
using IntPoint = typename Quantizer2D::IntPoint;
using CoordHash = typename Quantizer2D::CoordHash;

//------------------------------------------------------------------------------
// Helper: Create a square PLC from vertices
//------------------------------------------------------------------------------
PLC createSquarePLC() {
  PLC plc;

  // Square edges (4 edges forming a closed loop)
  plc.facets.resize(4);
  plc.facets[0] = {0, 1}; // bottom
  plc.facets[1] = {1, 2}; // right
  plc.facets[2] = {2, 3}; // top
  plc.facets[3] = {3, 0}; // left

  return plc;
}

//------------------------------------------------------------------------------
// Helper: Create square vertices with specified bounds
//------------------------------------------------------------------------------
vector<RealType> createSquareVertices(const RealPoint& lo = RealPoint(0.0, 0.0),
                                      const RealPoint& hi = RealPoint(1.0, 1.0)) {
  return {
    lo.x, lo.y,  // 0
    hi.x, lo.y,  // 1
    hi.x, hi.y,  // 2
    lo.x, hi.y   // 3
  };
}

//------------------------------------------------------------------------------
// Test: Basic construction and quantization
//------------------------------------------------------------------------------
void testBasicConstruction(const int tnum) {
  cout << "\n=== Test " << tnum << ": Basic Construction and Quantization ===" << endl;

  RealPoint xlo(0.0, 0.0);
  RealPoint xhi(1.0, 1.0);
  Quantizer2D Q(xlo, xhi);

  auto plc = createSquarePLC();
  auto vertices = createSquareVertices(xlo, xhi);

  QuantPLC2D qplc(plc, Q, vertices);

  // Check that all 4 vertices were quantized
  POLY_CHECK2(qplc.m_points.size() == 4,
              "Expected 4 vertices, got " << qplc.m_points.size());
  POLY_CHECK2(qplc.m_hashes.size() == 4,
              "Expected 4 hashes, got " << qplc.m_hashes.size());

  // Check that bounding box makes sense
  POLY_CHECK2(qplc.m_loBounds < qplc.m_hiBounds,
              "Invalid bounding box: lo >= hi");

  cout << "  Basic construction passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Quantization round-trip accuracy
//------------------------------------------------------------------------------
void testQuantizationAccuracy(const int tnum) {
  cout << "\n=== Test " << tnum << ": Quantization Accuracy ===" << endl;

  RealPoint xlo(-10.0, -10.0);
  RealPoint xhi(10.0, 10.0);
  Quantizer2D Q(xlo, xhi);

  vector<RealPoint> testPoints = {
    RealPoint(0.0, 0.0),
    RealPoint(5.0, 5.0),
    RealPoint(-5.0, -5.0),
    RealPoint(3.14, -2.71),
    RealPoint(-9.9, 9.9)
  };

  RealType maxError = Q.m_dx_o.x * 2.0;  // Allow 2x grid spacing

  for (const auto& p : testPoints) {
    IntPoint quantized = Q.quantize(p);
    RealPoint recovered = Q.dequantize(quantized);

    RealType dx = abs(recovered.x - p.x);
    RealType dy = abs(recovered.y - p.y);

    POLY_CHECK2(dx < maxError && dy < maxError,
                "Quantization error too large for point (" << p.x << ", " << p.y << ")");
  }

  cout << "  Quantization accuracy passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Deduplication via hashing (removeDegeneracies)
//------------------------------------------------------------------------------
void testDeduplication(const int tnum) {
  cout << "\n=== Test " << tnum << ": Hash-Based Deduplication ===" << endl;

  RealPoint xlo(0.0, 0.0);
  RealPoint xhi(1.0, 1.0);
  Quantizer2D Q(xlo, xhi);

  // Create vertices with duplicates (vertices 0 and 4 are the same)
  vector<RealType> vertices = {
    0.0, 0.0,  // 0
    1.0, 0.0,  // 1
    1.0, 1.0,  // 2
    0.0, 1.0,  // 3
    0.0, 0.0,  // 4 (duplicate of 0)
    1.0, 1.0   // 5 (duplicate of 2)
  };

  PLC plc;
  plc.facets.resize(2);
  plc.facets[0] = {0, 1};  // Uses 0
  plc.facets[1] = {4, 2};  // Uses 4 (duplicate of 0)

  QuantPLC2D qplc(plc, Q, vertices);

  // Before deduplication, should have 6 points
  POLY_CHECK2(qplc.m_points.size() == 6,
              "Expected 6 points before deduplication");

  // Remove duplicates
  qplc.removeDegeneracies();

  // After deduplication: 4 unique points (removed duplicates 4 and 5)
  POLY_CHECK2(qplc.m_points.size() == 4,
              "Expected 4 unique points after deduplication, got " << qplc.m_points.size());

  cout << "  Deduplication passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Reduce to only points on facets
//------------------------------------------------------------------------------
void testReduction(const int tnum) {
  cout << "\n=== Test " << tnum << ": Reduction (Unused Vertices) ===" << endl;

  RealPoint xlo(0.0, 0.0);
  RealPoint xhi(1.0, 1.0);
  Quantizer2D Q(xlo, xhi);

  // Create vertices with some unused
  vector<RealType> vertices = {
    0.0, 0.0,  // 0 - used
    1.0, 0.0,  // 1 - used
    1.0, 1.0,  // 2 - used
    0.5, 0.5,  // 3 - NOT used
    0.2, 0.7   // 4 - NOT used
  };

  PLC plc;
  plc.facets.resize(1);
  plc.facets[0] = {0, 1};  // Only uses vertices 0 and 1

  QuantPLC2D qplc(plc, Q, vertices);

  // Before reduction, should have 5 points
  POLY_CHECK2(qplc.m_points.size() == 5,
              "Expected 5 points before reduction");

  // Reduce: remove unused vertices
  qplc.reduce();

  // After reduction: 2 unique points (0, 1; removed unused 2, 3, 4)
  POLY_CHECK2(qplc.m_points.size() == 2,
              "Expected 2 unique points after reduction, got " << qplc.m_points.size());
  POLY_CHECK2(qplc.m_reduced, "Should be marked as reduced");

  cout << "  Reduction passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Convex hull computation (2D)
//------------------------------------------------------------------------------
void testConvexHull(const int tnum) {
  cout << "\n=== Test " << tnum << ": Convex Hull Computation ===" << endl;

  RealPoint xlo(-1.0, -1.0);
  RealPoint xhi(2.0, 2.0);
  Quantizer2D Q(xlo, xhi);

  // Create a set of points including some interior points
  vector<RealType> vertices = {
    0.0, 0.0,   // 0 (corner - should be in hull)
    1.0, 0.0,   // 1 (corner - should be in hull)
    1.0, 1.0,   // 2 (corner - should be in hull)
    0.0, 1.0,   // 3 (corner - should be in hull)
    0.5, 0.5,   // 4 (interior - should NOT be in hull)
    0.3, 0.7    // 5 (interior - should NOT be in hull)
  };

  PLC plc;  // Empty PLC, will compute hull
  QuantPLC2D qplc(plc, Q, vertices);

  qplc.makeConvex();

  // Square hull should have 4 edges
  POLY_CHECK2(qplc.facets.size() == 4,
              "Square hull should have 4 edges, got " << qplc.facets.size());

  // Each edge should have 2 vertices
  for (size_t i = 0; i < qplc.facets.size(); ++i) {
    POLY_CHECK2(qplc.facets[i].size() == 2,
                "Edge " << i << " should have 2 vertices, got " << qplc.facets[i].size());
  }

  // After reduction, should only have 4 points (the hull vertices)
  POLY_CHECK2(qplc.m_points.size() == 4,
              "Convex hull should have 4 vertices, got " << qplc.m_points.size());

  cout << "  Convex hull computation passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Edge ordering (facets form a closed loop)
//------------------------------------------------------------------------------
void testEdgeOrdering(const int tnum) {
  cout << "\n=== Test " << tnum << ": Edge Ordering ===" << endl;

  RealPoint xlo(0.0, 0.0);
  RealPoint xhi(1.0, 1.0);
  Quantizer2D Q(xlo, xhi);

  auto plc = createSquarePLC();
  auto vertices = createSquareVertices();

  QuantPLC2D qplc(plc, Q, vertices);
  qplc.reduce();  // This calls orderFacets

  // Check that edges form a closed loop
  for (size_t i = 0; i < qplc.facets.size(); ++i) {
    size_t next = (i + 1) % qplc.facets.size();
    POLY_CHECK2(qplc.facets[i][1] == qplc.facets[next][0],
                "Edge " << i << " should connect to edge " << next);
  }

  cout << "  Edge ordering passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Collinear vertex removal
//------------------------------------------------------------------------------
void testCollinearRemoval(const int tnum) {
  cout << "\n=== Test " << tnum << ": Collinear Vertex Removal ===" << endl;

  RealPoint xlo(0.0, 0.0);
  RealPoint xhi(10.0, 10.0);
  Quantizer2D Q(xlo, xhi);

  // Create a square with extra collinear points on edges
  vector<RealType> vertices = {
    0.0, 0.0,   // 0 - corner
    5.0, 0.0,   // 1 - collinear on bottom edge
    10.0, 0.0,  // 2 - corner
    10.0, 5.0,  // 3 - collinear on right edge
    10.0, 10.0, // 4 - corner
    5.0, 10.0,  // 5 - collinear on top edge
    0.0, 10.0,  // 6 - corner
    0.0, 5.0    // 7 - collinear on left edge
  };

  PLC plc;
  plc.facets.resize(8);
  plc.facets[0] = {0, 1};
  plc.facets[1] = {1, 2};
  plc.facets[2] = {2, 3};
  plc.facets[3] = {3, 4};
  plc.facets[4] = {4, 5};
  plc.facets[5] = {5, 6};
  plc.facets[6] = {6, 7};
  plc.facets[7] = {7, 0};

  QuantPLC2D qplc(plc, Q, vertices);

  // Before reduction, should have 8 edges
  POLY_CHECK2(qplc.facets.size() == 8,
              "Should have 8 edges before collinear removal");

  // Reduce calls orderFacets which removes collinear vertices
  qplc.reduce();

  // After collinear removal, should have 4 edges (square corners only)
  POLY_CHECK2(qplc.facets.size() == 4,
              "After removing collinear vertices, should have 4 edges, got " << qplc.facets.size());

  // Should have 4 unique vertices (the corners)
  POLY_CHECK2(qplc.m_points.size() == 4,
              "Should have 4 unique vertices, got " << qplc.m_points.size());

  cout << "  Collinear removal passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Point containment (within) - basic cases
//------------------------------------------------------------------------------
void testWithinBasic(const int tnum) {
  cout << "\n=== Test " << tnum << ": Point Containment (Basic) ===" << endl;

  RealPoint xlo(0.0, 0.0);
  RealPoint xhi(1.0, 1.0);
  Quantizer2D Q(xlo, xhi);

  auto plc = createSquarePLC();
  auto vertices = createSquareVertices(xlo, xhi);

  QuantPLC2D qplc(plc, Q, vertices);
  qplc.makeConvex();

  // Test points inside
  POLY_CHECK(qplc.within(RealPoint(0.5, 0.5)));  // Center
  POLY_CHECK(qplc.within(RealPoint(0.1, 0.1)));  // Near corner
  POLY_CHECK(qplc.within(RealPoint(0.9, 0.5)));  // Near edge

  // Test points outside
  POLY_CHECK(!qplc.within(RealPoint(-0.5, 0.5)));  // Left of square
  POLY_CHECK(!qplc.within(RealPoint(1.5, 0.5)));   // Right of square
  POLY_CHECK(!qplc.within(RealPoint(0.5, -0.5)));  // Below square
  POLY_CHECK(!qplc.within(RealPoint(0.5, 1.5)));   // Above square

  // Test points on boundary (should be inside)
  POLY_CHECK(qplc.within(RealPoint(0.0, 0.5)));  // Left edge
  POLY_CHECK(qplc.within(RealPoint(1.0, 0.5)));  // Right edge
  POLY_CHECK(qplc.within(RealPoint(0.5, 0.0)));  // Bottom edge
  POLY_CHECK(qplc.within(RealPoint(0.5, 1.0)));  // Top edge

  // Test vertices (should be inside)
  POLY_CHECK(qplc.within(RealPoint(0.0, 0.0)));
  POLY_CHECK(qplc.within(RealPoint(1.0, 1.0)));

  cout << "  Basic containment tests passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Point containment with holes
//------------------------------------------------------------------------------
void testWithinHoles(const int tnum) {
  cout << "\n=== Test " << tnum << ": Point Containment (With Holes) ===" << endl;

  RealPoint xlo(-10.0, -10.0);
  RealPoint xhi(10.0, 10.0);
  Quantizer2D Q(xlo, xhi);

  // Outer square: [-5, 5]^2 with inner hole: [-1, 1]^2
  vector<RealType> vertices = {
    // Outer square (0-3)
    -5.0, -5.0,  // 0
     5.0, -5.0,  // 1
     5.0,  5.0,  // 2
    -5.0,  5.0,  // 3
    // Inner square (4-7) - hole
    -1.0, -1.0,  // 4
     1.0, -1.0,  // 5
     1.0,  1.0,  // 6
    -1.0,  1.0   // 7
  };

  PLC plc;
  plc.facets.resize(4);
  plc.facets[0] = {0, 1};
  plc.facets[1] = {1, 2};
  plc.facets[2] = {2, 3};
  plc.facets[3] = {3, 0};

  plc.holes.resize(1);
  plc.holes[0].resize(4);
  plc.holes[0][0] = {4, 5};
  plc.holes[0][1] = {5, 6};
  plc.holes[0][2] = {6, 7};
  plc.holes[0][3] = {7, 4};

  QuantPLC2D qplc(plc, Q, vertices);

  // Inside outer, outside hole
  POLY_CHECK(qplc.within(RealPoint(-3.0, 0.0)));
  POLY_CHECK(qplc.within(RealPoint(3.0, 3.0)));

  // Inside hole (should be outside)
  POLY_CHECK(!qplc.within(RealPoint(0.0, 0.0)));
  POLY_CHECK(!qplc.within(RealPoint(0.5, 0.5)));

  // Outside outer boundary
  POLY_CHECK(!qplc.within(RealPoint(-6.0, 0.0)));
  POLY_CHECK(!qplc.within(RealPoint(0.0, 6.0)));

  // On hole boundary (should be inside by convention)
  POLY_CHECK(qplc.within(RealPoint(-1.0, 0.0)));
  POLY_CHECK(qplc.within(RealPoint(1.0, 0.0)));

  cout << "  Containment with holes passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Collinear function from GeomUtils
//------------------------------------------------------------------------------
void testCollinearFunction(const int tnum) {
  cout << "\n=== Test " << tnum << ": Collinear Function ===" << endl;

  // Test collinear points on a line segment
  IntPoint p1(0, 0);
  IntPoint p2(10, 0);
  IntPoint p_mid(5, 0);  // Midpoint, collinear
  IntPoint p_off(5, 1);  // Off the line

  POLY_CHECK(collinear(p1, p2, p_mid));   // Midpoint is collinear
  POLY_CHECK(!collinear(p1, p2, p_off));  // Off-line point is not collinear

  // Test point on line but outside segment
  IntPoint p_beyond(15, 0);  // Collinear but beyond p2
  POLY_CHECK(!collinear(p1, p2, p_beyond));  // Should be outside segment

  // Test diagonal line
  IntPoint d1(0, 0);
  IntPoint d2(10, 10);
  IntPoint d_mid(5, 5);
  IntPoint d_off(5, 6);

  POLY_CHECK(collinear(d1, d2, d_mid));
  POLY_CHECK(!collinear(d1, d2, d_off));

  cout << "  Collinear function passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Hash comparison utilities
//------------------------------------------------------------------------------
void testHashComparison(const int tnum) {
  cout << "\n=== Test " << tnum << ": Hash Comparison Utilities ===" << endl;

  RealPoint xlo(0.0, 0.0);
  RealPoint xhi(1.0, 1.0);
  Quantizer2D Q(xlo, xhi);

  auto plc1 = createSquarePLC();
  auto vertices1 = createSquareVertices();
  QuantPLC2D qplc1(plc1, Q, vertices1);

  // Create identical PLC
  auto plc2 = createSquarePLC();
  auto vertices2 = createSquareVertices();
  QuantPLC2D qplc2(plc2, Q, vertices2);

  // Should have same hashes (order-independent comparison)
  POLY_CHECK(QuantPLC2D::compareHashes(qplc1, qplc2));

  // Should have same facets (order-independent)
  POLY_CHECK(QuantPLC2D::compareFacets(qplc1, qplc2));

  cout << "  Hash comparison passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Stress test with many points
//------------------------------------------------------------------------------
void testStress(const int tnum) {
  cout << "\n=== Test " << tnum << ": Stress Test ===" << endl;

  RealPoint xlo(-100.0, -100.0);
  RealPoint xhi(100.0, 100.0);
  Quantizer2D Q(xlo, xhi);

  // Generate random point cloud
  const unsigned nPoints = 100;
  vector<RealType> vertices;
  vertices.reserve(nPoints * 2);

  for (unsigned i = 0; i < nPoints; ++i) {
    vertices.push_back((random01() - 0.5) * 200.0);  // x
    vertices.push_back((random01() - 0.5) * 200.0);  // y
  }

  PLC plc;  // Empty PLC
  QuantPLC2D qplc(plc, Q, vertices);

  // Compute convex hull
  qplc.makeConvex();

  // Hull should have at least 3 edges and at most n edges
  POLY_CHECK2(qplc.facets.size() >= 3,
              "Convex hull should have at least 3 edges");
  POLY_CHECK2(qplc.facets.size() <= nPoints,
              "Convex hull should have at most " << nPoints << " edges");

  // All edges should have 2 vertices
  for (size_t i = 0; i < qplc.facets.size(); ++i) {
    POLY_CHECK2(qplc.facets[i].size() == 2,
                "Edge " << i << " should have 2 vertices");
  }

  cout << "  Stress test passed! Hull has " << qplc.facets.size()
       << " edges from " << qplc.m_points.size() << " vertices" << endl;
}

} // anonymous namespace

//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Init(&argc, &argv);
#else
  POLY_CONTRACT_VAR(argc);
  POLY_CONTRACT_VAR(argv);
#endif

  srand(42);  // Deterministic randomness for reproducibility
  int tnum = 1;

  try {
    // Quantization tests
    testBasicConstruction(tnum++);
    testQuantizationAccuracy(tnum++);
    testDeduplication(tnum++);
    testReduction(tnum++);

    // Geometric tests
    testConvexHull(tnum++);
    testEdgeOrdering(tnum++);

    // Collinear removal tests
    testCollinearRemoval(tnum++);
    testCollinearFunction(tnum++);

    // Containment tests
    testWithinBasic(tnum++);
    testWithinHoles(tnum++);

    // Utility tests
    testHashComparison(tnum++);

    // Stress test
    testStress(tnum++);

    cout << "\n=== All QuantPLC2D tests passed! ===" << endl;

  } catch (const exception& e) {
    cerr << "\nTest failed with exception: " << e.what() << endl;
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
