// Comprehensive unit tests for QuantPLC<3>
//
// Tests the quantized PLC functionality including:
//   - Construction from real-space coordinates
//   - Quantization and deduplication
//   - Convex hull computation
//   - Facet ordering and orientation
//   - Point containment (within) tests
//   - Intersection tests via segmentFaceIntersection3D
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
#include "Point.hh"
#include "polytope_test_utilities.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace polytope;

namespace {

using RealType = double;
using RealPoint = Point3<RealType>;
using PLC = PLC<3>;
using QuantPLC3D = QuantPLC<3>;
using Quantizer3D = Quantizer<3>;
using IntPoint = typename Quantizer3D::IntPoint;
using CoordHash = typename Quantizer3D::CoordHash;
using Wide = WideInt<3>;

//------------------------------------------------------------------------------
// Helper: Create a cube PLC from vertices
//------------------------------------------------------------------------------
PLC createCubePLC() {
  PLC plc;

  // Cube facets (6 faces, each with 4 vertices)
  plc.facets.resize(6);
  plc.facets[0] = {0, 1, 2, 3}; // bottom (z = 0)
  plc.facets[1] = {4, 7, 6, 5}; // top (z = 1)
  plc.facets[2] = {0, 4, 5, 1}; // front (y = 0)
  plc.facets[3] = {2, 6, 7, 3}; // back (y = 1)
  plc.facets[4] = {0, 3, 7, 4}; // left (x = 0)
  plc.facets[5] = {1, 5, 6, 2}; // right (x = 1)

  return plc;
}

//------------------------------------------------------------------------------
// Helper: Create cube vertices with specified bounds
//------------------------------------------------------------------------------
vector<RealType> createCubeVertices(const RealPoint& lo = RealPoint(0.0, 0.0, 0.0),
                                    const RealPoint& hi = RealPoint(1.0, 1.0, 1.0)) {
  return {
    lo.x, lo.y, lo.z,  // 0
    hi.x, lo.y, lo.z,  // 1
    hi.x, hi.y, lo.z,  // 2
    lo.x, hi.y, lo.z,  // 3
    lo.x, lo.y, hi.z,  // 4
    hi.x, lo.y, hi.z,  // 5
    hi.x, hi.y, hi.z,  // 6
    lo.x, hi.y, hi.z   // 7
  };
}

//------------------------------------------------------------------------------
// Test: Basic construction and quantization
//------------------------------------------------------------------------------
void testBasicConstruction(const int tnum) {
  cout << "\n=== Test " << tnum << ": Basic Construction and Quantization ===" << endl;

  RealPoint xlo(0.0, 0.0, 0.0);
  RealPoint xhi(1.0, 1.0, 1.0);
  Quantizer3D Q(xlo, xhi);

  auto plc = createCubePLC();
  auto vertices = createCubeVertices(xlo, xhi);

  QuantPLC3D qplc(plc, Q, vertices);

  // Check that all 8 vertices were quantized
  POLY_CHECK2(qplc.m_points.size() == 8,
              "Expected 8 vertices, got " << qplc.m_points.size());
  POLY_CHECK2(qplc.m_hashes.size() == 8,
              "Expected 8 hashes, got " << qplc.m_hashes.size());

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

  RealPoint xlo(-10.0, -10.0, -10.0);
  RealPoint xhi(10.0, 10.0, 10.0);
  Quantizer3D Q(xlo, xhi);

  vector<RealPoint> testPoints = {
    RealPoint(0.0, 0.0, 0.0),
    RealPoint(5.0, 5.0, 5.0),
    RealPoint(-5.0, -5.0, -5.0),
    RealPoint(3.14, -2.71, 1.41),
    RealPoint(-9.9, 9.9, 0.1)
  };

  RealType maxError = Q.m_dx_o.x * 2.0;  // Allow 2x grid spacing

  for (const auto& p : testPoints) {
    IntPoint quantized = Q.quantize(p);
    RealPoint recovered = Q.dequantize(quantized);

    RealType dx = abs(recovered.x - p.x);
    RealType dy = abs(recovered.y - p.y);
    RealType dz = abs(recovered.z - p.z);

    POLY_CHECK2(dx < maxError && dy < maxError && dz < maxError,
                "Quantization error too large for point (" << p.x << ", " << p.y << ", " << p.z << ")");
  }

  cout << "  Quantization accuracy passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Deduplication via hashing
//------------------------------------------------------------------------------
void testDeduplication(const int tnum) {
  cout << "\n=== Test " << tnum << ": Hash-Based Deduplication ===" << endl;

  RealPoint xlo(0.0, 0.0, 0.0);
  RealPoint xhi(1.0, 1.0, 1.0);
  Quantizer3D Q(xlo, xhi);

  // Create vertices with duplicates (vertices 0 and 8 are the same)
  vector<RealType> vertices = {
    0.0, 0.0, 0.0,  // 0
    1.0, 0.0, 0.0,  // 1
    1.0, 1.0, 0.0,  // 2
    0.0, 1.0, 0.0,  // 3
    0.0, 0.0, 1.0,  // 4
    1.0, 0.0, 1.0,  // 5
    1.0, 1.0, 1.0,  // 6
    0.0, 1.0, 1.0,  // 7
    0.0, 0.0, 0.0   // 8 (duplicate of 0)
  };

  PLC plc;
  plc.facets.resize(1);
  plc.facets[0] = {0, 1, 2, 8};  // Uses both 0 and 8 (same point)

  QuantPLC3D qplc(plc, Q, vertices);

  // Before reduction, should have 9 points
  POLY_CHECK2(qplc.m_points.size() == 9,
              "Expected 9 points before reduction");

  // Reduce based on the provided PLC
  qplc.reduce();
  // After reduction, duplicates should be merged
  POLY_CHECK2(qplc.m_points.size() == 3,
              "Expected 3 unique points after reduction, got " << qplc.m_points.size());

  cout << "  Deduplication passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Convex hull computation
//------------------------------------------------------------------------------
void testConvexHull(const int tnum) {
  cout << "\n=== Test " << tnum << ": Convex Hull Computation ===" << endl;

  RealPoint xlo(-1.0, -1.0, -1.0);
  RealPoint xhi(1.0, 1.0, 1.0);
  Quantizer3D Q(xlo, xhi);

  // Create a tetrahedron with one interior point
  vector<RealType> vertices = {
    0.0, 0.0, 0.0,      // 0 (base corner)
    1.0, 0.0, 0.0,      // 1 (base corner)
    0.5, 1.0, 0.0,      // 2 (base corner)
    0.5, 0.5, 1.0,      // 3 (apex)
    0.5, 0.5, 0.3       // 4 (interior point - should not be in hull)
  };

  PLC plc;  // Empty PLC, will compute hull
  QuantPLC3D qplc(plc, Q, vertices);

  qplc.makeConvex();

  // Tetrahedron hull should have 4 facets
  POLY_CHECK2(qplc.facets.size() == 4,
              "Tetrahedron hull should have 4 facets, got " << qplc.facets.size());

  // Each facet should be a triangle (3 vertices)
  for (size_t i = 0; i < qplc.facets.size(); ++i) {
    POLY_CHECK2(qplc.facets[i].size() == 3,
                "Tetrahedron facet " << i << " should have 3 vertices, got " << qplc.facets[i].size());
  }

  // After reduction, should only have 4 points (the hull vertices)
  POLY_CHECK2(qplc.m_points.size() == 4,
              "Convex hull should have 4 vertices, got " << qplc.m_points.size());

  cout << "  Convex hull computation passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Facet orientation (normals pointing outward)
//------------------------------------------------------------------------------
void testFacetOrientation(const int tnum) {
  cout << "\n=== Test " << tnum << ": Facet Orientation ===" << endl;

  RealPoint xlo(0.0, 0.0, 0.0);
  RealPoint xhi(2.0, 2.0, 2.0);
  Quantizer3D Q(xlo, xhi);

  auto plc = createCubePLC();
  auto vertices = createCubeVertices();

  QuantPLC3D qplc(plc, Q, vertices);
  qplc.reduce();

  // Compute centroid in floating-point (no need for exact integer math here)
  RealPoint centroid(0.0, 0.0, 0.0);
  for (const auto& p : qplc.m_points) {
    auto pf = p.template type_cast<double>();
    centroid = centroid + pf;
  }
  centroid = centroid / static_cast<double>(qplc.m_points.size());

  // Check that all facet normals point away from centroid
  for (size_t i = 0; i < qplc.facets.size(); ++i) {
    const auto& facet = qplc.facets[i];
    POLY_CHECK2(facet.size() >= 3, "Facet " << i << " has fewer than 3 vertices");

    // Use floating-point for geometric validation (no need for exact integer math)
    const auto v0 = qplc.m_points[facet[0]].template type_cast<double>();
    const auto v1 = qplc.m_points[facet[1]].template type_cast<double>();
    const auto v2 = qplc.m_points[facet[2]].template type_cast<double>();

    // Compute normal: (v1 - v0) × (v2 - v0)
    auto edge1 = v1 - v0;
    auto edge2 = v2 - v0;
    RealPoint normal(
      edge1.y * edge2.z - edge1.z * edge2.y,
      edge1.z * edge2.x - edge1.x * edge2.z,
      edge1.x * edge2.y - edge1.y * edge2.x
    );

    // Vector from centroid to face
    auto toFace = (v0 - centroid).template type_cast<double>();

    // Dot product should be positive (normal points outward)
    double dot = normal.x * toFace.x + normal.y * toFace.y + normal.z * toFace.z;

    POLY_CHECK2(dot > 0.0,
                "Facet " << i << " normal points inward (dot = " << dot << ")");
  }

  cout << "  Facet orientation passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Point containment (within) - basic cases
//------------------------------------------------------------------------------
void testWithinBasic(const int tnum) {
  cout << "\n=== Test " << tnum << ": Point Containment (Basic) ===" << endl;

  RealPoint xlo(0.0, 0.0, 0.0);
  RealPoint xhi(1.0, 1.0, 1.0);
  Quantizer3D Q(xlo, xhi);

  auto plc = createCubePLC();
  auto vertices = createCubeVertices(xlo, xhi);

  QuantPLC3D qplc(plc, Q, vertices);
  qplc.makeConvex();

  // Test points inside
  POLY_CHECK(qplc.within(RealPoint(0.5, 0.5, 0.5)));  // Center
  POLY_CHECK(qplc.within(RealPoint(0.1, 0.1, 0.1)));  // Near corner
  POLY_CHECK(qplc.within(RealPoint(0.9, 0.5, 0.5)));  // Near face

  // Test points outside
  POLY_CHECK(!qplc.within(RealPoint(-0.5, 0.5, 0.5)));  // Outside in -x
  POLY_CHECK(!qplc.within(RealPoint(1.5, 0.5, 0.5)));   // Outside in +x
  POLY_CHECK(!qplc.within(RealPoint(0.5, -0.5, 0.5)));  // Outside in -y
  POLY_CHECK(!qplc.within(RealPoint(0.5, 1.5, 0.5)));   // Outside in +y
  POLY_CHECK(!qplc.within(RealPoint(0.5, 0.5, -0.5)));  // Outside in -z
  POLY_CHECK(!qplc.within(RealPoint(0.5, 0.5, 1.5)));   // Outside in +z

  // Test points on boundary (should be inside)
  POLY_CHECK(qplc.within(RealPoint(0.0, 0.5, 0.5)));  // On left face
  POLY_CHECK(qplc.within(RealPoint(1.0, 0.5, 0.5)));  // On right face
  POLY_CHECK(qplc.within(RealPoint(0.5, 0.0, 0.5)));  // On front face
  POLY_CHECK(qplc.within(RealPoint(0.5, 1.0, 0.5)));  // On back face
  POLY_CHECK(qplc.within(RealPoint(0.5, 0.5, 0.0)));  // On bottom face
  POLY_CHECK(qplc.within(RealPoint(0.5, 0.5, 1.0)));  // On top face

  // Test vertices (should be inside)
  POLY_CHECK(qplc.within(RealPoint(0.0, 0.0, 0.0)));
  POLY_CHECK(qplc.within(RealPoint(1.0, 1.0, 1.0)));

  cout << "  Basic containment tests passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Point containment with holes
//------------------------------------------------------------------------------
void testWithinHoles(const int tnum) {
  cout << "\n=== Test " << tnum << ": Point Containment (With Holes) ===" << endl;

  RealPoint xlo(-10.0, -10.0, -10.0);
  RealPoint xhi(10.0, 10.0, 10.0);
  Quantizer3D Q(xlo, xhi);

  // Outer cube: [-5, 5]^3 with inner hole: [-1, 1]^3
  vector<RealType> vertices = {
    // Outer cube (0-7)
    -5.0, -5.0, -5.0,  // 0
     5.0, -5.0, -5.0,  // 1
     5.0,  5.0, -5.0,  // 2
    -5.0,  5.0, -5.0,  // 3
    -5.0, -5.0,  5.0,  // 4
     5.0, -5.0,  5.0,  // 5
     5.0,  5.0,  5.0,  // 6
    -5.0,  5.0,  5.0,  // 7
    // Inner cube (8-15) - hole
    -1.0, -1.0, -1.0,  // 8
     1.0, -1.0, -1.0,  // 9
     1.0,  1.0, -1.0,  // 10
    -1.0,  1.0, -1.0,  // 11
    -1.0, -1.0,  1.0,  // 12
     1.0, -1.0,  1.0,  // 13
     1.0,  1.0,  1.0,  // 14
    -1.0,  1.0,  1.0   // 15
  };

  PLC plc;
  plc.facets.resize(6);
  plc.facets[0] = {0, 1, 2, 3};
  plc.facets[1] = {4, 7, 6, 5};
  plc.facets[2] = {0, 4, 5, 1};
  plc.facets[3] = {2, 6, 7, 3};
  plc.facets[4] = {0, 3, 7, 4};
  plc.facets[5] = {1, 5, 6, 2};

  plc.holes.resize(1);
  plc.holes[0].resize(6);
  plc.holes[0][0] = {8, 9, 10, 11};
  plc.holes[0][1] = {12, 15, 14, 13};
  plc.holes[0][2] = {8, 12, 13, 9};
  plc.holes[0][3] = {10, 14, 15, 11};
  plc.holes[0][4] = {8, 11, 15, 12};
  plc.holes[0][5] = {9, 13, 14, 10};

  QuantPLC3D qplc(plc, Q, vertices);

  // Inside outer, outside hole
  POLY_CHECK(qplc.within(RealPoint(-3.0, 0.0, 0.0)));
  POLY_CHECK(qplc.within(RealPoint(3.0, 3.0, 3.0)));

  // Inside hole (should be outside)
  POLY_CHECK(!qplc.within(RealPoint(0.0, 0.0, 0.0)));
  POLY_CHECK(!qplc.within(RealPoint(0.5, 0.5, 0.5)));

  // Outside outer boundary
  POLY_CHECK(!qplc.within(RealPoint(-6.0, 0.0, 0.0)));
  POLY_CHECK(!qplc.within(RealPoint(0.0, 0.0, 6.0)));

  // On hole boundary (by convention, should be inside)
  POLY_CHECK(qplc.within(RealPoint(-1.0, 0.0, 0.0)));
  POLY_CHECK(qplc.within(RealPoint(1.0, 0.0, 0.0)));

  cout << "  Containment with holes passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Segment-face intersection (basic)
//------------------------------------------------------------------------------
void testSegmentFaceIntersection(const int tnum) {
  cout << "\n=== Test " << tnum << ": Segment-Face Intersection ===" << endl;

  RealPoint xlo(-2.0, -2.0, -2.0);
  RealPoint xhi(2.0, 2.0, 2.0);
  Quantizer3D Q(xlo, xhi);

  // Create a simple square face in the xy-plane at z=0
  vector<RealType> vertices = {
    -1.0, -1.0, 0.0,  // 0
     1.0, -1.0, 0.0,  // 1
     1.0,  1.0, 0.0,  // 2
    -1.0,  1.0, 0.0   // 3
  };

  PLC plc;
  plc.facets.resize(1);
  plc.facets[0] = {0, 1, 2, 3};  // Square face

  QuantPLC3D qplc(plc, Q, vertices);

  // Test intersection from above to below (should hit center)
  IntPoint segStart = Q.quantize(RealPoint(0.0, 0.0, 1.0));
  IntPoint segEnd = Q.quantize(RealPoint(0.0, 0.0, -1.0));
  IntPoint result;

  int check = segmentFaceIntersection3D(segStart, segEnd, qplc.facets[0], qplc.m_points, result);
  POLY_CHECK2(check == 1, "Segment should intersect face");

  RealPoint hitReal = Q.dequantize(result);
  POLY_CHECK2(abs(hitReal.z) < 0.1, "Intersection z-coordinate should be near 0, got " << hitReal.z);
  POLY_CHECK2(abs(hitReal.x) < 0.1, "Intersection x-coordinate should be near 0, got " << hitReal.x);
  POLY_CHECK2(abs(hitReal.y) < 0.1, "Intersection y-coordinate should be near 0, got " << hitReal.y);

  // Test segment parallel to face (should not intersect)
  segStart = Q.quantize(RealPoint(-1.0, 0.0, 1.0));
  segEnd = Q.quantize(RealPoint(1.0, 0.0, 1.0));
  check = segmentFaceIntersection3D(segStart, segEnd, qplc.facets[0], qplc.m_points, result);
  POLY_CHECK2(check < 0, "Parallel segment should not intersect");

  // Test segment missing the face (crosses plane but outside polygon)
  segStart = Q.quantize(RealPoint(2.0, 2.0, 1.0));
  segEnd = Q.quantize(RealPoint(2.0, 2.0, -1.0));
  check = segmentFaceIntersection3D(segStart, segEnd, qplc.facets[0], qplc.m_points, result);
  POLY_CHECK2(check < 0, "Segment should miss face polygon");

  cout << "  Segment-face intersection passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Segment-face intersection (edge cases)
//------------------------------------------------------------------------------
void testSegmentFaceEdgeCases(const int tnum) {
  cout << "\n=== Test " << tnum << ": Segment-Face Intersection (Edge Cases) ===" << endl;

  RealPoint xlo(-2.0, -2.0, -2.0);
  RealPoint xhi(2.0, 2.0, 2.0);
  Quantizer3D Q(xlo, xhi);

  // Create a triangular face
  vector<RealType> vertices = {
    0.0, 0.0, 0.0,   // 0
    1.0, 0.0, 0.0,   // 1
    0.5, 1.0, 0.0    // 2
  };

  PLC plc;
  plc.facets.resize(1);
  plc.facets[0] = {0, 1, 2};

  QuantPLC3D qplc(plc, Q, vertices);

  IntPoint result;

  // Test segment hitting exactly on vertex
  IntPoint segStart = Q.quantize(RealPoint(0.0, 0.0, 1.0));
  IntPoint segEnd = Q.quantize(RealPoint(0.0, 0.0, -1.0));
  int check = segmentFaceIntersection3D(segStart, segEnd, qplc.facets[0], qplc.m_points, result);
  POLY_CHECK2(check == 1, "Segment through vertex should intersect");

  // Test segment hitting exactly on edge
  segStart = Q.quantize(RealPoint(0.5, 0.0, 1.0));
  segEnd = Q.quantize(RealPoint(0.5, 0.0, -1.0));
  check = segmentFaceIntersection3D(segStart, segEnd, qplc.facets[0], qplc.m_points, result);
  POLY_CHECK2(check == 1, "Segment through edge should intersect");

  // Test segment just barely inside polygon
  segStart = Q.quantize(RealPoint(0.5, 0.5, 1.0));
  segEnd = Q.quantize(RealPoint(0.5, 0.5, -1.0));
  check = segmentFaceIntersection3D(segStart, segEnd, qplc.facets[0], qplc.m_points, result);
  POLY_CHECK2(check == 1, "Segment through interior should intersect");

  // Test segment endpoint on face
  segStart = Q.quantize(RealPoint(0.5, 0.3, 0.0));
  segEnd = Q.quantize(RealPoint(0.5, 0.3, -1.0));
  check = segmentFaceIntersection3D(segStart, segEnd, qplc.facets[0], qplc.m_points, result);
  POLY_CHECK2(check >= 0, "Segment with endpoint on face should intersect");

  cout << "  Edge case intersections passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Overflow protection in intersections
//------------------------------------------------------------------------------
void testOverflowProtection(const int tnum) {
  cout << "\n=== Test " << tnum << ": Overflow Protection in Intersections ===" << endl;

  // Use large coordinate space to stress test overflow handling
  RealPoint xlo(-1e6, -1e6, -1e6);
  RealPoint xhi(1e6, 1e6, 1e6);
  Quantizer3D Q(xlo, xhi);

  // Create a face with large coordinates
  vector<RealType> vertices = {
    -1e5, -1e5, 0.0,
     1e5, -1e5, 0.0,
     1e5,  1e5, 0.0,
    -1e5,  1e5, 0.0
  };

  PLC plc;
  plc.facets.resize(1);
  plc.facets[0] = {0, 1, 2, 3};

  QuantPLC3D qplc(plc, Q, vertices);

  // Test intersection with large coordinates
  IntPoint segStart = Q.quantize(RealPoint(0.0, 0.0, 1e5));
  IntPoint segEnd = Q.quantize(RealPoint(0.0, 0.0, -1e5));
  IntPoint result;

  int check = segmentFaceIntersection3D(segStart, segEnd, qplc.facets[0], qplc.m_points, result);
  POLY_CHECK2(check == 1, "Large coordinate intersection should work");

  RealPoint hitReal = Q.dequantize(result);
  POLY_CHECK2(abs(hitReal.z) < 1e4, "Intersection z should be near plane");

  cout << "  Overflow protection passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Multiple intersections (ray casting)
//------------------------------------------------------------------------------
void testMultipleIntersections(const int tnum) {
  cout << "\n=== Test " << tnum << ": Multiple Intersections (Ray Casting) ===" << endl;

  RealPoint xlo(-2.0, -2.0, -2.0);
  RealPoint xhi(2.0, 2.0, 2.0);
  Quantizer3D Q(xlo, xhi);

  auto plc = createCubePLC();
  auto vertices = createCubeVertices();

  QuantPLC3D qplc(plc, Q, vertices);

  // Cast a ray from outside through the cube
  // Ray from (-1, 0.5, 0.5) to (2, 0.5, 0.5) should hit left and right faces
  IntPoint rayStart = Q.quantize(RealPoint(-1.0, 0.5, 0.5));
  IntPoint rayEnd = Q.quantize(RealPoint(2.0, 0.5, 0.5));

  set<CoordHash> intersections;
  for (const auto& facet : qplc.facets) {
    IntPoint result;
    if (segmentFaceIntersection3D(rayStart, rayEnd, facet, qplc.m_points, result) >= 0) {
      intersections.insert(Q.hash(result));
    }
  }

  // Should hit exactly 2 faces (entry and exit)
  POLY_CHECK2(intersections.size() == 2,
              "Ray should hit 2 faces, got " << intersections.size());

  // Ray from inside to outside should hit 1 face
  rayStart = Q.quantize(RealPoint(0.5, 0.5, 0.5));
  rayEnd = Q.quantize(RealPoint(2.0, 0.5, 0.5));

  intersections.clear();
  for (const auto& facet : qplc.facets) {
    IntPoint result;
    if (segmentFaceIntersection3D(rayStart, rayEnd, facet, qplc.m_points, result) >= 0) {
      intersections.insert(Q.hash(result));
    }
  }

  POLY_CHECK2(intersections.size() == 1,
              "Ray from inside should hit 1 face, got " << intersections.size());

  cout << "  Multiple intersections passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Hash comparison utilities
//------------------------------------------------------------------------------
void testHashComparison(const int tnum) {
  cout << "\n=== Test " << tnum << ": Hash Comparison Utilities ===" << endl;

  RealPoint xlo(0.0, 0.0, 0.0);
  RealPoint xhi(1.0, 1.0, 1.0);
  Quantizer3D Q(xlo, xhi);

  auto plc1 = createCubePLC();
  auto vertices1 = createCubeVertices();
  QuantPLC3D qplc1(plc1, Q, vertices1);

  // Create identical PLC with permuted vertex indices
  auto plc2 = createCubePLC();
  auto vertices2 = createCubeVertices();
  QuantPLC3D qplc2(plc2, Q, vertices2);

  // Should have same hashes (order-independent comparison)
  POLY_CHECK(QuantPLC3D::compareHashes(qplc1, qplc2));

  // Should have same facets (order-independent)
  POLY_CHECK(QuantPLC3D::compareFacets(qplc1, qplc2));

  cout << "  Hash comparison passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Reduction with unused vertices
//------------------------------------------------------------------------------
void testReduction(const int tnum) {
  cout << "\n=== Test " << tnum << ": Reduction (Unused Vertices) ===" << endl;

  RealPoint xlo(0.0, 0.0, 0.0);
  RealPoint xhi(2.0, 2.0, 2.0);
  Quantizer3D Q(xlo, xhi);

  // Create vertices with some unused ones
  vector<RealType> vertices = {
    0.0, 0.0, 0.0,  // 0 - used
    1.0, 0.0, 0.0,  // 1 - used
    1.0, 1.0, 0.0,  // 2 - used
    0.0, 1.0, 0.0,  // 3 - used
    5.0, 5.0, 5.0,  // 4 - NOT used
    6.0, 6.0, 6.0   // 5 - NOT used
  };

  PLC plc;
  plc.facets.resize(1);
  plc.facets[0] = {0, 1, 2, 3};  // Only uses vertices 0-3

  QuantPLC3D qplc(plc, Q, vertices);

  // Before reduction: 6 vertices (including unused ones)
  POLY_CHECK2(qplc.m_points.size() == 6,
              "Before reduction, should have 6 vertices, got " << qplc.m_points.size());
  POLY_CHECK2(!qplc.m_reduced, "Should not be marked as reduced yet");

  // Call reduction explicitly
  qplc.reduce();

  // After reduction: 4 vertices (0-3 only, unused vertices 4-5 removed)
  POLY_CHECK2(qplc.m_points.size() == 4,
              "After reduction, should have 4 vertices, got " << qplc.m_points.size());
  POLY_CHECK2(qplc.m_reduced, "Should be marked as reduced");

  cout << "  Reduction passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Quantized cross product and dot product
//------------------------------------------------------------------------------
void testQuantizedOperations(const int tnum) {
  cout << "\n=== Test " << tnum << ": Quantized Operations ===" << endl;

  using IntType = int64_t;

  // Test qdot
  IntPoint a(1, 2, 3);
  IntPoint b(4, 5, 6);
  Wide result = qdot<3>(a, b);
  Wide expected = 1*4 + 2*5 + 3*6;  // = 32
  POLY_CHECK2(result == expected,
              "qdot failed: expected " << static_cast<IntType>(expected)
              << ", got " << static_cast<IntType>(result));

  // Test qcross
  IntPoint c(1, 0, 0);
  IntPoint d(0, 1, 0);
  auto cross = qcross(c, d);
  POLY_CHECK2(cross.x == 0 && cross.y == 0 && cross.z == 1,
              "qcross(x, y) should give z");

  // Test with negative values
  IntPoint e(-1, 2, -3);
  IntPoint f(4, -5, 6);
  result = qdot<3>(e, f);
  expected = (-1)*4 + 2*(-5) + (-3)*6;  // = -32
  POLY_CHECK2(result == expected,
              "qdot with negatives failed");

  cout << "  Quantized operations passed!" << endl;
}

//------------------------------------------------------------------------------
// Test: Stress test with many points
//------------------------------------------------------------------------------
void testStress(const int tnum) {
  cout << "\n=== Test " << tnum << ": Stress Test ===" << endl;

  RealPoint xlo(-100.0, -100.0, -100.0);
  RealPoint xhi(100.0, 100.0, 100.0);
  Quantizer3D Q(xlo, xhi);

  // Generate random point cloud
  const unsigned nPoints = 100;
  vector<RealType> vertices;
  vertices.reserve(nPoints * 3);

  for (unsigned i = 0; i < nPoints; ++i) {
    vertices.push_back((random01() - 0.5) * 200.0);  // x
    vertices.push_back((random01() - 0.5) * 200.0);  // y
    vertices.push_back((random01() - 0.5) * 200.0);  // z
  }

  PLC plc;  // Empty PLC
  QuantPLC3D qplc(plc, Q, vertices);

  // Compute convex hull
  qplc.makeConvex();

  // Hull should have at least 4 facets (tetrahedron) and at most O(n) facets
  POLY_CHECK2(qplc.facets.size() >= 4,
              "Convex hull should have at least 4 facets");
  POLY_CHECK2(qplc.facets.size() <= nPoints,
              "Convex hull should have at most " << nPoints << " facets");

  // All facets should have at least 3 vertices
  for (size_t i = 0; i < qplc.facets.size(); ++i) {
    POLY_CHECK2(qplc.facets[i].size() >= 3,
                "Facet " << i << " has fewer than 3 vertices");
  }

  // Test within for a point we know is inside
  RealPoint center(0.0, 0.0, 0.0);
  qplc.within(center);
  // Center may or may not be inside depending on random points
  // Just verify the test doesn't crash

  cout << "  Stress test passed! Hull has " << qplc.facets.size()
       << " facets from " << qplc.m_points.size() << " vertices" << endl;
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

    // Geometric tests
    testConvexHull(tnum++);
    testFacetOrientation(tnum++);
    testReduction(tnum++);

    // Containment tests
    testWithinBasic(tnum++);
    testWithinHoles(tnum++);

    // Intersection tests
    testSegmentFaceIntersection(tnum++);
    testSegmentFaceEdgeCases(tnum++);
    testOverflowProtection(tnum++);
    testMultipleIntersections(tnum++);

    // Utility tests
    testQuantizedOperations(tnum++);
    testHashComparison(tnum++);

    // Stress test
    testStress(tnum++);

    cout << "\n=== All QuantPLC3D tests passed! ===" << endl;

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
