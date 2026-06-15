// Comprehensive unit tests for 3D convex intersection routines
//
// Tests the quantized convex intersection functionality in Intersections.hh:
//   - Basic cube-cube intersection tests
//   - Non-intersecting polyhedra (separated by various axes)
//   - Edge cases: touching faces, shared vertices, shared edges
//   - Separating axis theorem validation
//   - Point containment vs edge-edge separation
//   - Degenerate cases
//   - Stress tests with complex polyhedra

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

#include "polytope.hh"
#include "QuantPLC.hh"
#include "Quantizer.hh"
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

//------------------------------------------------------------------------------
// Helper: Create a cube PLC with 6 facets
//------------------------------------------------------------------------------
PLC createCubePLC() {
  PLC plc;

  // Cube facets (6 faces, each with 4 vertices)
  plc.facets.resize(6);
  plc.facets[0] = {0, 1, 2, 3}; // bottom (z = lo)
  plc.facets[1] = {4, 7, 6, 5}; // top (z = hi)
  plc.facets[2] = {0, 4, 5, 1}; // front (y = lo)
  plc.facets[3] = {2, 6, 7, 3}; // back (y = hi)
  plc.facets[4] = {0, 3, 7, 4}; // left (x = lo)
  plc.facets[5] = {1, 5, 6, 2}; // right (x = hi)

  return plc;
}

//------------------------------------------------------------------------------
// Helper: Create cube vertices with specified bounds
//------------------------------------------------------------------------------
std::vector<RealType> createCubeVertices(const RealPoint& lo = RealPoint(0.0, 0.0, 0.0),
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
// Helper: Create a tetrahedron PLC
//------------------------------------------------------------------------------
PLC createTetrahedronPLC() {
  PLC plc;

  // Tetrahedron has 4 facets
  plc.facets.resize(4);
  plc.facets[0] = {0, 1, 2};  // base
  plc.facets[1] = {0, 3, 1};  // side
  plc.facets[2] = {1, 3, 2};  // side
  plc.facets[3] = {2, 3, 0};  // side

  return plc;
}

//------------------------------------------------------------------------------
// Helper: Create tetrahedron vertices centered at origin with given size
//------------------------------------------------------------------------------
std::vector<RealType> createTetrahedronVertices(const RealPoint& center = RealPoint(0.0, 0.0, 0.0),
                                                RealType size = 1.0) {
  return {
    center.x - size, center.y - size, center.z - size,  // 0
    center.x + size, center.y - size, center.z - size,  // 1
    center.x,        center.y + size, center.z - size,  // 2
    center.x,        center.y,        center.z + size   // 3 (apex)
  };
}

//------------------------------------------------------------------------------
// Test 1: Two overlapping cubes
//------------------------------------------------------------------------------
void testOverlappingCubes(const int tnum) {
  std::cout << "\n=== Test " << tnum << ": Overlapping Cubes ===" << std::endl;

  RealPoint xlo(-2.0, -2.0, -2.0);
  RealPoint xhi(2.0, 2.0, 2.0);
  Quantizer3D Q(xlo, xhi);

  // Cube A: [0, 1]^3
  auto plcA = createCubePLC();
  auto verticesA = createCubeVertices(RealPoint(0.0, 0.0, 0.0), RealPoint(1.0, 1.0, 1.0));
  QuantPLC3D qplcA(plcA, Q, verticesA);

  // Cube B: [0.5, 1.5]^3 (overlaps with A)
  auto plcB = createCubePLC();
  auto verticesB = createCubeVertices(RealPoint(0.5, 0.5, 0.5), RealPoint(1.5, 1.5, 1.5));
  QuantPLC3D qplcB(plcB, Q, verticesB);

  // Make convex (computes normals)
  qplcA.makeConvex();
  qplcB.makeConvex();

  // Test intersection
  bool intersects = QuantPLC3D::convexPLCIntersection(qplcA, qplcB);

  POLY_CHECK2(intersects, "Overlapping cubes should intersect");

  std::cout << "  Overlapping cubes test passed!" << std::endl;
}

//------------------------------------------------------------------------------
// Test 2: Two separated cubes (no intersection)
//------------------------------------------------------------------------------
void testSeparatedCubes(const int tnum) {
  std::cout << "\n=== Test " << tnum << ": Separated Cubes ===" << std::endl;

  RealPoint xlo(-5.0, -5.0, -5.0);
  RealPoint xhi(5.0, 5.0, 5.0);
  Quantizer3D Q(xlo, xhi);

  // Cube A: [0, 1]^3
  auto plcA = createCubePLC();
  auto verticesA = createCubeVertices(RealPoint(0.0, 0.0, 0.0), RealPoint(1.0, 1.0, 1.0));
  QuantPLC3D qplcA(plcA, Q, verticesA);

  // Cube B: [2, 3]^3 (separated from A in all dimensions)
  auto plcB = createCubePLC();
  auto verticesB = createCubeVertices(RealPoint(2.0, 2.0, 2.0), RealPoint(3.0, 3.0, 3.0));
  QuantPLC3D qplcB(plcB, Q, verticesB);

  // Make convex (computes normals)
  qplcA.makeConvex();
  qplcB.makeConvex();

  // Test intersection
  bool intersects = QuantPLC3D::convexPLCIntersection(qplcA, qplcB);

  POLY_CHECK2(!intersects, "Separated cubes should not intersect");

  std::cout << "  Separated cubes test passed!" << std::endl;
}

//------------------------------------------------------------------------------
// Test 3: Cubes touching at a face
//------------------------------------------------------------------------------
void testTouchingFaces(const int tnum) {
  std::cout << "\n=== Test " << tnum << ": Touching Faces ===" << std::endl;

  RealPoint xlo(-2.0, -2.0, -2.0);
  RealPoint xhi(2.0, 2.0, 2.0);
  Quantizer3D Q(xlo, xhi);

  // Cube A: [0, 1]^3
  auto plcA = createCubePLC();
  auto verticesA = createCubeVertices(RealPoint(0.0, 0.0, 0.0), RealPoint(1.0, 1.0, 1.0));
  QuantPLC3D qplcA(plcA, Q, verticesA);

  // Cube B: [1, 2]^3 (shares the x=1 face with A)
  auto plcB = createCubePLC();
  auto verticesB = createCubeVertices(RealPoint(1.0, 0.0, 0.0), RealPoint(2.0, 1.0, 1.0));
  QuantPLC3D qplcB(plcB, Q, verticesB);

  // Make convex (computes normals)
  qplcA.makeConvex();
  qplcB.makeConvex();

  // Test intersection (touching counts as intersection)
  bool intersects = QuantPLC3D::convexPLCIntersection(qplcA, qplcB);

  POLY_CHECK2(intersects, "Cubes touching at a face should intersect");

  std::cout << "  Touching faces test passed!" << std::endl;
}

//------------------------------------------------------------------------------
// Test 4: Cubes touching at an edge
//------------------------------------------------------------------------------
void testTouchingEdges(const int tnum) {
  std::cout << "\n=== Test " << tnum << ": Touching Edges ===" << std::endl;

  RealPoint xlo(-2.0, -2.0, -2.0);
  RealPoint xhi(2.0, 2.0, 2.0);
  Quantizer3D Q(xlo, xhi);

  // Cube A: [0, 1]^3
  auto plcA = createCubePLC();
  auto verticesA = createCubeVertices(RealPoint(0.0, 0.0, 0.0), RealPoint(1.0, 1.0, 1.0));
  QuantPLC3D qplcA(plcA, Q, verticesA);

  // Cube B: [1, 2] x [1, 2] x [0, 1] (shares edge along z-axis at x=1, y=1)
  auto plcB = createCubePLC();
  auto verticesB = createCubeVertices(RealPoint(1.0, 1.0, 0.0), RealPoint(2.0, 2.0, 1.0));
  QuantPLC3D qplcB(plcB, Q, verticesB);

  // Make convex (computes normals)
  qplcA.makeConvex();
  qplcB.makeConvex();

  // Test intersection (touching at edge counts as intersection)
  bool intersects = QuantPLC3D::convexPLCIntersection(qplcA, qplcB);

  POLY_CHECK2(intersects, "Cubes touching at an edge should intersect");

  std::cout << "  Touching edges test passed!" << std::endl;
}

//------------------------------------------------------------------------------
// Test 5: Cubes touching at a vertex
//------------------------------------------------------------------------------
void testTouchingVertex(const int tnum) {
  std::cout << "\n=== Test " << tnum << ": Touching Vertex ===" << std::endl;

  RealPoint xlo(-2.0, -2.0, -2.0);
  RealPoint xhi(2.0, 2.0, 2.0);
  Quantizer3D Q(xlo, xhi);

  // Cube A: [0, 1]^3
  auto plcA = createCubePLC();
  auto verticesA = createCubeVertices(RealPoint(0.0, 0.0, 0.0), RealPoint(1.0, 1.0, 1.0));
  QuantPLC3D qplcA(plcA, Q, verticesA);

  // Cube B: [1, 2]^3 (shares only corner vertex (1,1,1))
  auto plcB = createCubePLC();
  auto verticesB = createCubeVertices(RealPoint(1.0, 1.0, 1.0), RealPoint(2.0, 2.0, 2.0));
  QuantPLC3D qplcB(plcB, Q, verticesB);

  // Make convex (computes normals)
  qplcA.makeConvex();
  qplcB.makeConvex();

  // Test intersection (touching at vertex counts as intersection)
  bool intersects = QuantPLC3D::convexPLCIntersection(qplcA, qplcB);

  POLY_CHECK2(intersects, "Cubes touching at a vertex should intersect");

  std::cout << "  Touching vertex test passed!" << std::endl;
}

//------------------------------------------------------------------------------
// Test 6: One polyhedron contains the other
//------------------------------------------------------------------------------
void testContainment(const int tnum) {
  std::cout << "\n=== Test " << tnum << ": One Contains the Other ===" << std::endl;

  RealPoint xlo(-5.0, -5.0, -5.0);
  RealPoint xhi(5.0, 5.0, 5.0);
  Quantizer3D Q(xlo, xhi);

  // Large cube: [-2, 2]^3
  auto plcA = createCubePLC();
  auto verticesA = createCubeVertices(RealPoint(-2.0, -2.0, -2.0), RealPoint(2.0, 2.0, 2.0));
  QuantPLC3D qplcA(plcA, Q, verticesA);

  // Small cube inside: [-1, 1]^3
  auto plcB = createCubePLC();
  auto verticesB = createCubeVertices(RealPoint(-1.0, -1.0, -1.0), RealPoint(1.0, 1.0, 1.0));
  QuantPLC3D qplcB(plcB, Q, verticesB);

  // Make convex (computes normals)
  qplcA.makeConvex();
  qplcB.makeConvex();

  // Test intersection (containment is intersection)
  bool intersects = QuantPLC3D::convexPLCIntersection(qplcA, qplcB);

  POLY_CHECK2(intersects, "Large cube containing small cube should intersect");

  std::cout << "  Containment test passed!" << std::endl;
}

//------------------------------------------------------------------------------
// Test 7: Tetrahedra intersection
//------------------------------------------------------------------------------
void testTetrahedraIntersect(const int tnum) {
  std::cout << "\n=== Test " << tnum << ": Overlapping Tetrahedra ===" << std::endl;

  RealPoint xlo(-3.0, -3.0, -3.0);
  RealPoint xhi(3.0, 3.0, 3.0);
  Quantizer3D Q(xlo, xhi);

  // Tetrahedron A centered at origin
  auto plcA = createTetrahedronPLC();
  auto verticesA = createTetrahedronVertices(RealPoint(0.0, 0.0, 0.0), 1.0);
  QuantPLC3D qplcA(plcA, Q, verticesA);

  // Tetrahedron B offset slightly (overlaps with A)
  auto plcB = createTetrahedronPLC();
  auto verticesB = createTetrahedronVertices(RealPoint(0.5, 0.5, 0.5), 1.0);
  QuantPLC3D qplcB(plcB, Q, verticesB);

  // Make convex (computes normals)
  qplcA.makeConvex();
  qplcB.makeConvex();

  // Test intersection
  bool intersects = QuantPLC3D::convexPLCIntersection(qplcA, qplcB);

  POLY_CHECK2(intersects, "Overlapping tetrahedra should intersect");

  std::cout << "  Overlapping tetrahedra test passed!" << std::endl;
}

//------------------------------------------------------------------------------
// Test 8: Tetrahedra non-intersecting
//------------------------------------------------------------------------------
void testTetrahedraSeparated(const int tnum) {
  std::cout << "\n=== Test " << tnum << ": Separated Tetrahedra ===" << std::endl;

  RealPoint xlo(-5.0, -5.0, -5.0);
  RealPoint xhi(5.0, 5.0, 5.0);
  Quantizer3D Q(xlo, xhi);

  // Tetrahedron A centered at origin
  auto plcA = createTetrahedronPLC();
  auto verticesA = createTetrahedronVertices(RealPoint(0.0, 0.0, 0.0), 1.0);
  QuantPLC3D qplcA(plcA, Q, verticesA);

  // Tetrahedron B far away
  auto plcB = createTetrahedronPLC();
  auto verticesB = createTetrahedronVertices(RealPoint(3.0, 3.0, 3.0), 1.0);
  QuantPLC3D qplcB(plcB, Q, verticesB);

  // Make convex (computes normals)
  qplcA.makeConvex();
  qplcB.makeConvex();

  // Test intersection
  bool intersects = QuantPLC3D::convexPLCIntersection(qplcA, qplcB);

  POLY_CHECK2(!intersects, "Separated tetrahedra should not intersect");

  std::cout << "  Separated tetrahedra test passed!" << std::endl;
}

//------------------------------------------------------------------------------
// Test 9: Cube and tetrahedron intersection
//------------------------------------------------------------------------------
void testCubeTetrahedronIntersect(const int tnum) {
  std::cout << "\n=== Test " << tnum << ": Cube-Tetrahedron Intersection ===" << std::endl;

  RealPoint xlo(-3.0, -3.0, -3.0);
  RealPoint xhi(3.0, 3.0, 3.0);
  Quantizer3D Q(xlo, xhi);

  // Cube: [0, 1]^3
  auto plcA = createCubePLC();
  auto verticesA = createCubeVertices(RealPoint(0.0, 0.0, 0.0), RealPoint(1.0, 1.0, 1.0));
  QuantPLC3D qplcA(plcA, Q, verticesA);

  // Tetrahedron overlapping cube
  auto plcB = createTetrahedronPLC();
  auto verticesB = createTetrahedronVertices(RealPoint(0.5, 0.5, 0.5), 0.8);
  QuantPLC3D qplcB(plcB, Q, verticesB);

  // Make convex (computes normals)
  qplcA.makeConvex();
  qplcB.makeConvex();

  // Test intersection
  bool intersects = QuantPLC3D::convexPLCIntersection(qplcA, qplcB);

  POLY_CHECK2(intersects, "Cube and tetrahedron should intersect");

  std::cout << "  Cube-tetrahedron intersection test passed!" << std::endl;
}

//------------------------------------------------------------------------------
// Test 10: Cube and tetrahedron separated
//------------------------------------------------------------------------------
void testCubeTetrahedronSeparated(const int tnum) {
  std::cout << "\n=== Test " << tnum << ": Cube-Tetrahedron Separated ===" << std::endl;

  RealPoint xlo(-5.0, -5.0, -5.0);
  RealPoint xhi(5.0, 5.0, 5.0);
  Quantizer3D Q(xlo, xhi);

  // Cube: [0, 1]^3
  auto plcA = createCubePLC();
  auto verticesA = createCubeVertices(RealPoint(0.0, 0.0, 0.0), RealPoint(1.0, 1.0, 1.0));
  QuantPLC3D qplcA(plcA, Q, verticesA);

  // Tetrahedron far away
  auto plcB = createTetrahedronPLC();
  auto verticesB = createTetrahedronVertices(RealPoint(3.0, 3.0, 3.0), 0.8);
  QuantPLC3D qplcB(plcB, Q, verticesB);

  // Make convex (computes normals)
  qplcA.makeConvex();
  qplcB.makeConvex();

  // Test intersection
  bool intersects = QuantPLC3D::convexPLCIntersection(qplcA, qplcB);

  POLY_CHECK2(!intersects, "Separated cube and tetrahedron should not intersect");

  std::cout << "  Cube-tetrahedron separation test passed!" << std::endl;
}

//------------------------------------------------------------------------------
// Test 11: Rotated cubes (separating axis test)
//------------------------------------------------------------------------------
void testRotatedCubes(const int tnum) {
  std::cout << "\n=== Test " << tnum << ": Rotated Cubes ===" << std::endl;

  RealPoint xlo(-3.0, -3.0, -3.0);
  RealPoint xhi(3.0, 3.0, 3.0);
  Quantizer3D Q(xlo, xhi);

  // Cube A: axis-aligned [0, 1]^3
  auto plcA = createCubePLC();
  auto verticesA = createCubeVertices(RealPoint(0.0, 0.0, 0.0), RealPoint(1.0, 1.0, 1.0));
  QuantPLC3D qplcA(plcA, Q, verticesA);

  // Cube B: rotated 45 degrees around z-axis, centered at (1.5, 0.5, 0.5)
  // This creates a diamond shape in xy-plane
  RealType s = 0.707;  // sin(45°) ≈ cos(45°)
  std::vector<RealType> verticesB = {
    1.5 - s, 0.5,     0.0,   // 0
    1.5,     0.5 - s, 0.0,   // 1
    1.5 + s, 0.5,     0.0,   // 2
    1.5,     0.5 + s, 0.0,   // 3
    1.5 - s, 0.5,     1.0,   // 4
    1.5,     0.5 - s, 1.0,   // 5
    1.5 + s, 0.5,     1.0,   // 6
    1.5,     0.5 + s, 1.0    // 7
  };

  auto plcB = createCubePLC();
  QuantPLC3D qplcB(plcB, Q, verticesB);

  // Make convex (computes normals)
  qplcA.makeConvex();
  qplcB.makeConvex();

  // These cubes overlap slightly at the boundary
  bool intersects = QuantPLC3D::convexPLCIntersection(qplcA, qplcB);

  POLY_CHECK2(intersects, "Rotated cubes with overlap should intersect");

  std::cout << "  Rotated cubes test passed!" << std::endl;
}

//------------------------------------------------------------------------------
// Test 11: Rotated cubes (separating axis test)
//------------------------------------------------------------------------------
void testShiftedRotatedCubes(const int tnum) {
  std::cout << "\n=== Test " << tnum << ": Shifted Rotated Cubes ===" << std::endl;

  RealPoint xlo(-3.0, -3.0, -3.0);
  RealPoint xhi(4.0, 4.0, 4.0);
  Quantizer3D Q(xlo, xhi);

  // Cube A: axis-aligned [0, 1]^3
  auto plcA = createCubePLC();
  auto verticesA = createCubeVertices(RealPoint(0.0, 0.0, 0.0), RealPoint(1.0, 1.0, 1.0));
  QuantPLC3D qplcA(plcA, Q, verticesA);

  // Cube B: rotated 45 degrees around z-axis, centered at (1.5, 0.5, 0.5)
  // This creates a diamond shape in xy-plane
  RealType s = 0.707;  // sin(45°) ≈ cos(45°)
  std::vector<RealType> verticesB = {
    2.5 - s, 0.5,     0.0,   // 0
    2.5,     0.5 - s, 0.0,   // 1
    2.5 + s, 0.5,     0.0,   // 2
    2.5,     0.5 + s, 0.0,   // 3
    2.5 - s, 0.5,     1.0,   // 4
    2.5,     0.5 - s, 1.0,   // 5
    2.5 + s, 0.5,     1.0,   // 6
    2.5,     0.5 + s, 1.0    // 7
  };

  auto plcB = createCubePLC();
  QuantPLC3D qplcB(plcB, Q, verticesB);

  // Make convex (computes normals)
  qplcA.makeConvex();
  qplcB.makeConvex();

  // These cubes overlap slightly at the boundary
  bool intersects = QuantPLC3D::convexPLCIntersection(qplcA, qplcB);

  POLY_CHECK2(!intersects, "Shifted rotated cubes without overlap should not intersect");

  std::cout << "  Shifted rotated cubes test passed!" << std::endl;
}

//------------------------------------------------------------------------------
// Test 12: Elongated polyhedra (edge-edge separation)
//------------------------------------------------------------------------------
void testElongatedPolyhedra(const int tnum) {
  std::cout << "\n=== Test " << tnum << ": Elongated Polyhedra ===" << std::endl;

  RealPoint xlo(-5.0, -5.0, -5.0);
  RealPoint xhi(5.0, 5.0, 5.0);
  Quantizer3D Q(xlo, xhi);

  // Elongated box A along x-axis: [0, 3] x [0, 0.5] x [0, 0.5]
  auto plcA = createCubePLC();
  auto verticesA = createCubeVertices(RealPoint(0.0, 0.0, 0.0), RealPoint(3.0, 0.5, 0.5));
  QuantPLC3D qplcA(plcA, Q, verticesA);

  // Elongated box B along y-axis: [1, 1.5] x [0, 3] x [0, 0.5]
  auto plcB = createCubePLC();
  auto verticesB = createCubeVertices(RealPoint(1.0, 0.0, 0.0), RealPoint(1.5, 3.0, 0.5));
  QuantPLC3D qplcB(plcB, Q, verticesB);

  // Make convex (computes normals)
  qplcA.makeConvex();
  qplcB.makeConvex();

  // These should intersect (both pass through region around (1.25, 0.25, 0.25))
  bool intersects = QuantPLC3D::convexPLCIntersection(qplcA, qplcB);

  POLY_CHECK2(intersects, "Elongated polyhedra crossing should intersect");

  std::cout << "  Elongated polyhedra test passed!" << std::endl;
}

//------------------------------------------------------------------------------
// Test 13: Near-miss separation (stress test for separating axis)
//------------------------------------------------------------------------------
void testNearMissSeparation(const int tnum) {
  std::cout << "\n=== Test " << tnum << ": Near-Miss Separation ===" << std::endl;

  RealPoint xlo(-2.0, -2.0, -2.0);
  RealPoint xhi(2.0, 2.0, 2.0);
  Quantizer3D Q(xlo, xhi);

  // Cube A: [0, 1]^3
  auto plcA = createCubePLC();
  auto verticesA = createCubeVertices(RealPoint(0.0, 0.0, 0.0), RealPoint(1.0, 1.0, 1.0));
  QuantPLC3D qplcA(plcA, Q, verticesA);

  // Cube B: [1.01, 2.01]^3 (very close but separated by 0.01)
  auto plcB = createCubePLC();
  auto verticesB = createCubeVertices(RealPoint(1.01, 0.0, 0.0), RealPoint(2.01, 1.0, 1.0));
  QuantPLC3D qplcB(plcB, Q, verticesB);

  // Make convex (computes normals)
  qplcA.makeConvex();
  qplcB.makeConvex();

  // Should NOT intersect (separated by small gap)
  bool intersects = QuantPLC3D::convexPLCIntersection(qplcA, qplcB);

  POLY_CHECK2(!intersects, "Cubes separated by small gap should not intersect");

  std::cout << "  Near-miss separation test passed!" << std::endl;
}

//------------------------------------------------------------------------------
// Test 14: Symmetry test (A intersects B iff B intersects A)
//------------------------------------------------------------------------------
void testSymmetry(const int tnum) {
  std::cout << "\n=== Test " << tnum << ": Intersection Symmetry ===" << std::endl;

  RealPoint xlo(-3.0, -3.0, -3.0);
  RealPoint xhi(3.0, 3.0, 3.0);
  Quantizer3D Q(xlo, xhi);

  // Create various pairs of polyhedra
  std::vector<std::pair<QuantPLC3D, QuantPLC3D>> testPairs;

  // Pair 1: Overlapping cubes
  {
    auto plcA = createCubePLC();
    auto verticesA = createCubeVertices(RealPoint(0.0, 0.0, 0.0), RealPoint(1.0, 1.0, 1.0));
    QuantPLC3D qplcA(plcA, Q, verticesA);

    auto plcB = createCubePLC();
    auto verticesB = createCubeVertices(RealPoint(0.5, 0.5, 0.5), RealPoint(1.5, 1.5, 1.5));
    QuantPLC3D qplcB(plcB, Q, verticesB);

    testPairs.push_back({qplcA, qplcB});
  }

  // Pair 2: Separated cubes
  {
    auto plcA = createCubePLC();
    auto verticesA = createCubeVertices(RealPoint(0.0, 0.0, 0.0), RealPoint(1.0, 1.0, 1.0));
    QuantPLC3D qplcA(plcA, Q, verticesA);

    auto plcB = createCubePLC();
    auto verticesB = createCubeVertices(RealPoint(2.0, 2.0, 2.0), RealPoint(3.0, 3.0, 3.0));
    QuantPLC3D qplcB(plcB, Q, verticesB);

    testPairs.push_back({qplcA, qplcB});
  }

  // Pair 3: Cube and tetrahedron
  {
    auto plcA = createCubePLC();
    auto verticesA = createCubeVertices(RealPoint(0.0, 0.0, 0.0), RealPoint(1.0, 1.0, 1.0));
    QuantPLC3D qplcA(plcA, Q, verticesA);

    auto plcB = createTetrahedronPLC();
    auto verticesB = createTetrahedronVertices(RealPoint(0.5, 0.5, 0.5), 0.8);
    QuantPLC3D qplcB(plcB, Q, verticesB);

    testPairs.push_back({qplcA, qplcB});
  }

  // Make all test pairs convex
  for (auto& pair : testPairs) {
    pair.first.makeConvex();
    pair.second.makeConvex();
  }

  // Test symmetry for all pairs
  for (size_t i = 0; i < testPairs.size(); ++i) {
    const auto& qplcA = testPairs[i].first;
    const auto& qplcB = testPairs[i].second;

    bool AB = QuantPLC3D::convexPLCIntersection(qplcA, qplcB);
    bool BA = QuantPLC3D::convexPLCIntersection(qplcB, qplcA);

    POLY_CHECK2(AB == BA,
                "Symmetry failed for pair " << i << ": A∩B = " << AB << ", B∩A = " << BA);
  }

  std::cout << "  Intersection symmetry test passed!" << std::endl;
}

//------------------------------------------------------------------------------
// Test 15: Degenerate cases
//------------------------------------------------------------------------------
void testDegenerateCases(const int tnum) {
  std::cout << "\n=== Test " << tnum << ": Degenerate Cases ===" << std::endl;

  RealPoint xlo(-2.0, -2.0, -2.0);
  RealPoint xhi(2.0, 2.0, 2.0);
  Quantizer3D Q(xlo, xhi);

  // Test with same polyhedron (should always intersect with itself)
  auto plcA = createCubePLC();
  auto verticesA = createCubeVertices(RealPoint(0.0, 0.0, 0.0), RealPoint(1.0, 1.0, 1.0));
  QuantPLC3D qplcA(plcA, Q, verticesA);
  qplcA.makeConvex();

  bool intersects = QuantPLC3D::convexPLCIntersection(qplcA, qplcA);

  POLY_CHECK2(intersects, "Polyhedron should intersect with itself");

  std::cout << "  Degenerate cases test passed!" << std::endl;
}

//------------------------------------------------------------------------------
// Test 16: Stress test with complex convex hulls
//------------------------------------------------------------------------------
void testComplexHulls(const int tnum) {
  std::cout << "\n=== Test " << tnum << ": Complex Convex Hulls ===" << std::endl;

  RealPoint xlo(-10.0, -10.0, -10.0);
  RealPoint xhi(10.0, 10.0, 10.0);
  Quantizer3D Q(xlo, xhi);

  // Generate random point cloud A and compute its convex hull
  const unsigned nPointsA = 20;
  std::vector<RealType> verticesA;
  verticesA.reserve(nPointsA * 3);
  for (unsigned i = 0; i < nPointsA; ++i) {
    verticesA.push_back((random01() - 0.5) * 4.0);  // x in [-2, 2]
    verticesA.push_back((random01() - 0.5) * 4.0);  // y in [-2, 2]
    verticesA.push_back((random01() - 0.5) * 4.0);  // z in [-2, 2]
  }

  PLC plcA;  // Empty PLC
  QuantPLC3D qplcA(plcA, Q, verticesA);
  qplcA.makeConvex();

  // Generate random point cloud B (offset to ensure overlap)
  const unsigned nPointsB = 20;
  std::vector<RealType> verticesB;
  verticesB.reserve(nPointsB * 3);
  for (unsigned i = 0; i < nPointsB; ++i) {
    verticesB.push_back((random01() - 0.5) * 4.0 + 1.0);  // x in [-1, 3]
    verticesB.push_back((random01() - 0.5) * 4.0);        // y in [-2, 2]
    verticesB.push_back((random01() - 0.5) * 4.0);        // z in [-2, 2]
  }

  PLC plcB;  // Empty PLC
  QuantPLC3D qplcB(plcB, Q, verticesB);
  qplcB.makeConvex();

  // Test intersection (should likely intersect due to overlap)
  bool intersects = QuantPLC3D::convexPLCIntersection(qplcA, qplcB);

  // Just verify it doesn't crash - result depends on random points
  std::cout << "  Complex hulls test completed (intersects = " << intersects << ")" << std::endl;
  std::cout << "  Hull A: " << qplcA.m_points.size() << " vertices, "
            << qplcA.facets.size() << " facets" << std::endl;
  std::cout << "  Hull B: " << qplcB.m_points.size() << " vertices, "
            << qplcB.facets.size() << " facets" << std::endl;
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
    // Basic intersection tests
    testOverlappingCubes(tnum++);
    testSeparatedCubes(tnum++);

    // Boundary touching tests
    testTouchingFaces(tnum++);
    testTouchingEdges(tnum++);
    testTouchingVertex(tnum++);

    // Containment test
    testContainment(tnum++);

    // Tetrahedra tests
    testTetrahedraIntersect(tnum++);
    testTetrahedraSeparated(tnum++);

    // Mixed polyhedra tests
    testCubeTetrahedronIntersect(tnum++);
    testCubeTetrahedronSeparated(tnum++);

    // Separating axis theorem tests
    testRotatedCubes(tnum++);
    testShiftedRotatedCubes(tnum++);
    testElongatedPolyhedra(tnum++);
    testNearMissSeparation(tnum++);

    // Property tests
    testSymmetry(tnum++);
    testDegenerateCases(tnum++);

    // Stress test
    testComplexHulls(tnum++);

    std::cout << "\n=== All 3D convex intersection tests passed! ===" << std::endl;

  } catch (const std::exception& e) {
    std::cerr << "\nTest failed with exception: " << e.what() << std::endl;
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
