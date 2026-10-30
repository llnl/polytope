//------------------------------------------------------------------------------
// Test for 3D polyhedron-plane clipping routines
//
// Tests the layered clipping architecture:
//   Layer 1: Point classification (classifyPointByPlane)
//   Layer 2: Edge clipping (clipEdgeByPlane)
//   Layer 3: Face clipping (clipFaceByPlane)
//   Layer 4: Polyhedron clipping (clipPolyhedronByPlane)
//------------------------------------------------------------------------------

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

#include "Intersections.hh"
#include "Point.hh"

using namespace polytope;

//------------------------------------------------------------------------------
// Helper: Print a point
//------------------------------------------------------------------------------
template<typename CoordType>
void printPoint(const Point3<CoordType>& p, const std::string& label = "") {
  if (!label.empty()) std::cout << label << ": ";
  std::cout << "(" << p.x << ", " << p.y << ", " << p.z << ")" << std::endl;
}

//------------------------------------------------------------------------------
// Helper: Create a unit cube centered at origin
// Returns: (vertices, faces as vertex indices)
//------------------------------------------------------------------------------
template<typename CoordType>
std::pair<std::vector<Point3<CoordType>>, std::vector<std::vector<unsigned>>>
createUnitCube() {
  std::vector<Point3<CoordType>> vertices = {
    Point3<CoordType>(-1, -1, -1),  // 0
    Point3<CoordType>( 1, -1, -1),  // 1
    Point3<CoordType>( 1,  1, -1),  // 2
    Point3<CoordType>(-1,  1, -1),  // 3
    Point3<CoordType>(-1, -1,  1),  // 4
    Point3<CoordType>( 1, -1,  1),  // 5
    Point3<CoordType>( 1,  1,  1),  // 6
    Point3<CoordType>(-1,  1,  1)   // 7
  };

  // Faces with outward-pointing normals (right-hand rule)
  std::vector<std::vector<unsigned>> faces = {
    {0, 1, 2, 3},  // bottom (-z)
    {4, 7, 6, 5},  // top (+z)
    {0, 4, 5, 1},  // front (-y)
    {2, 6, 7, 3},  // back (+y)
    {0, 3, 7, 4},  // left (-x)
    {1, 5, 6, 2}   // right (+x)
  };

  return {vertices, faces};
}

//------------------------------------------------------------------------------
// Test 1: Layer 1 - Point classification
//------------------------------------------------------------------------------
void testPointClassification() {
  std::cout << "\n=== Test 1: Point Classification ===" << std::endl;

  using CoordType = int;
  Point3<CoordType> planePoint(0, 0, 0);
  Point3<CoordType> planeNormal(0, 0, 1);  // z-axis, pointing up

  Point3<CoordType> above(0, 0, 5);
  Point3<CoordType> below(0, 0, -5);
  Point3<CoordType> on(0, 0, 0);
  Point3<CoordType> on2(3, 4, 0);

  int classAbove = classifyPointByPlane(above, planePoint, planeNormal);
  int classBelow = classifyPointByPlane(below, planePoint, planeNormal);
  int classOn = classifyPointByPlane(on, planePoint, planeNormal);
  int classOn2 = classifyPointByPlane(on2, planePoint, planeNormal);

  std::cout << "Point (0,0,5) relative to z=0 plane: " << classAbove
            << " (expect +1)" << std::endl;
  std::cout << "Point (0,0,-5) relative to z=0 plane: " << classBelow
            << " (expect -1)" << std::endl;
  std::cout << "Point (0,0,0) relative to z=0 plane: " << classOn
            << " (expect 0)" << std::endl;
  std::cout << "Point (3,4,0) relative to z=0 plane: " << classOn2
            << " (expect 0)" << std::endl;

  assert(classAbove == 1);
  assert(classBelow == -1);
  assert(classOn == 0);
  assert(classOn2 == 0);

  std::cout << "PASS" << std::endl;
}

//------------------------------------------------------------------------------
// Test 2: Layer 2 - Edge clipping
//------------------------------------------------------------------------------
void testEdgeClipping() {
  std::cout << "\n=== Test 2: Edge Clipping ===" << std::endl;

  using CoordType = int;
  Point3<CoordType> planePoint(0, 0, 0);
  Point3<CoordType> planeNormal(0, 0, 1);  // z=0 plane

  // Edge crossing plane
  Point3<CoordType> a(-5, 0, -10);
  Point3<CoordType> b(5, 0, 10);
  Point3<CoordType> result;

  int status = clipEdgeByPlane(a, b, planePoint, planeNormal, result);

  std::cout << "Edge from (-5,0,-10) to (5,0,10) clips at: ";
  printPoint(result, "");
  std::cout << "Status: " << status << " (expect 1)" << std::endl;
  std::cout << "Expected: (0, 0, 0)" << std::endl;

  assert(status == 1);
  assert(result.x == 0 && result.y == 0 && result.z == 0);

  // Edge entirely above plane
  Point3<CoordType> c(0, 0, 5);
  Point3<CoordType> d(0, 0, 10);
  status = clipEdgeByPlane(c, d, planePoint, planeNormal, result);

  std::cout << "Edge from (0,0,5) to (0,0,10) status: " << status
            << " (expect -1, no crossing)" << std::endl;
  assert(status == -1);

  // Edge with start on plane
  Point3<CoordType> e(0, 0, 0);
  Point3<CoordType> f(0, 0, 10);
  status = clipEdgeByPlane(e, f, planePoint, planeNormal, result);

  std::cout << "Edge starting on plane: ";
  printPoint(result, "");
  std::cout << "Status: " << status << " (expect 0, start on plane)" << std::endl;

  assert(status == 0);
  assert(result.x == 0 && result.y == 0 && result.z == 0);

  std::cout << "PASS" << std::endl;
}

//------------------------------------------------------------------------------
// Test 3: Layer 3 - Face clipping
//------------------------------------------------------------------------------
void testFaceClipping() {
  std::cout << "\n=== Test 3: Face Clipping ===" << std::endl;

  using CoordType = int;
  Point3<CoordType> planePoint(0, 0, 0);
  Point3<CoordType> planeNormal(0, 0, 1);  // z=0 plane, keep z>0

  // Square face straddling the plane
  std::vector<Point3<CoordType>> faceVerts = {
    Point3<CoordType>(-5, -5, -5),
    Point3<CoordType>( 5, -5, -5),
    Point3<CoordType>( 5,  5,  5),
    Point3<CoordType>(-5,  5,  5)
  };

  auto result = clipFaceByPlane(faceVerts, planePoint, planeNormal, true);

  std::cout << "Original face: 4 vertices" << std::endl;
  std::cout << "Clipped face: " << result.vertices.size() << " vertices (expect 4)" << std::endl;
  std::cout << "Fully clipped: " << result.fullyClipped << " (expect 0)" << std::endl;
  std::cout << "Fully retained: " << result.fullyRetained << " (expect 0)" << std::endl;

  std::cout << "Clipped vertices:" << std::endl;
  for (size_t i = 0; i < result.vertices.size(); ++i) {
    std::cout << "  " << i << ": ";
    printPoint(result.vertices[i], "");
  }

  assert(!result.fullyClipped);
  assert(!result.fullyRetained);
  assert(result.vertices.size() == 4);  // Should create a clipped quadrilateral

  // Face entirely above plane (should be retained)
  std::vector<Point3<CoordType>> faceAbove = {
    Point3<CoordType>(-5, -5, 5),
    Point3<CoordType>( 5, -5, 5),
    Point3<CoordType>( 5,  5, 5),
    Point3<CoordType>(-5,  5, 5)
  };

  result = clipFaceByPlane(faceAbove, planePoint, planeNormal, true);

  std::cout << "\nFace entirely above plane:" << std::endl;
  std::cout << "Fully retained: " << result.fullyRetained << " (expect 1)" << std::endl;
  assert(result.fullyRetained);

  // Face entirely below plane (should be clipped away)
  std::vector<Point3<CoordType>> faceBelow = {
    Point3<CoordType>(-5, -5, -5),
    Point3<CoordType>( 5, -5, -5),
    Point3<CoordType>( 5,  5, -5),
    Point3<CoordType>(-5,  5, -5)
  };

  result = clipFaceByPlane(faceBelow, planePoint, planeNormal, true);

  std::cout << "\nFace entirely below plane:" << std::endl;
  std::cout << "Fully clipped: " << result.fullyClipped << " (expect 1)" << std::endl;
  assert(result.fullyClipped);

  std::cout << "PASS" << std::endl;
}

//------------------------------------------------------------------------------
// Test 4: Layer 4 - Polyhedron clipping
//------------------------------------------------------------------------------
void testPolyhedronClipping() {
  std::cout << "\n=== Test 4: Polyhedron Clipping ===" << std::endl;

  using CoordType = int;

  // Create a unit cube
  auto [vertices, faces] = createUnitCube<CoordType>();

  std::cout << "Original cube: " << vertices.size() << " vertices, "
            << faces.size() << " faces" << std::endl;

  // Clip by z=0 plane, keeping z>0 (should cut cube in half)
  Point3<CoordType> planePoint(0, 0, 0);
  Point3<CoordType> planeNormal(0, 0, 1);

  auto result = clipPolyhedronByPlane(vertices, faces, planePoint, planeNormal);

  std::cout << "Clipped result: " << result.vertices.size() << " vertices, "
            << result.faces.size() << " faces" << std::endl;
  std::cout << "Fully clipped: " << result.fullyClipped << " (expect 0)" << std::endl;
  std::cout << "Fully retained: " << result.fullyRetained << " (expect 0)" << std::endl;

  assert(!result.fullyClipped);
  assert(!result.fullyRetained);
  assert(result.vertices.size() >= 4);  // Should have at least the 4 top vertices
  assert(result.faces.size() >= 5);     // Top + 4 sides + cap

  std::cout << "\nClipped vertices:" << std::endl;
  for (size_t i = 0; i < result.vertices.size(); ++i) {
    std::cout << "  " << i << ": ";
    printPoint(result.vertices[i], "");
  }

  std::cout << "\nClipped faces:" << std::endl;
  for (size_t i = 0; i < result.faces.size(); ++i) {
    std::cout << "  Face " << i << ": [";
    for (size_t j = 0; j < result.faces[i].size(); ++j) {
      std::cout << result.faces[i][j];
      if (j < result.faces[i].size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
  }

  // Test: cube entirely above plane (should be retained)
  std::vector<Point3<CoordType>> verticesAbove;
  for (const auto& v : vertices) {
    verticesAbove.push_back(Point3<CoordType>(v.x, v.y, v.z + 10));
  }

  result = clipPolyhedronByPlane(verticesAbove, faces, planePoint, planeNormal);
  std::cout << "\nCube entirely above plane:" << std::endl;
  std::cout << "Fully retained: " << result.fullyRetained << " (expect 1)" << std::endl;
  assert(result.fullyRetained);

  // Test: cube entirely below plane (should be clipped away)
  std::vector<Point3<CoordType>> verticesBelow;
  for (const auto& v : vertices) {
    verticesBelow.push_back(Point3<CoordType>(v.x, v.y, v.z - 10));
  }

  result = clipPolyhedronByPlane(verticesBelow, faces, planePoint, planeNormal);
  std::cout << "\nCube entirely below plane:" << std::endl;
  std::cout << "Fully clipped: " << result.fullyClipped << " (expect 1)" << std::endl;
  assert(result.fullyClipped);

  std::cout << "PASS" << std::endl;
}

//------------------------------------------------------------------------------
// Test 5: Iterative clipping (simulate PLC clipping)
//------------------------------------------------------------------------------
void testIterativeClipping() {
  std::cout << "\n=== Test 5: Iterative Clipping ===" << std::endl;

  using CoordType = int;

  // Create a unit cube
  auto [vertices, faces] = createUnitCube<CoordType>();

  std::cout << "Starting with cube: " << vertices.size() << " vertices" << std::endl;

  // Define a sequence of clipping planes (simulating a complex PLC boundary)
  std::vector<std::pair<Point3<CoordType>, Point3<CoordType>>> planes = {
    {{0, 0, 0}, {0, 0, 1}},   // z = 0, keep z > 0
    {{0, 0, 0}, {1, 0, 0}},   // x = 0, keep x > 0
    {{0, 0, 0}, {0, 1, 0}}    // y = 0, keep y > 0
  };

  // Apply each plane iteratively
  auto currentVertices = vertices;
  auto currentFaces = faces;

  for (size_t i = 0; i < planes.size(); ++i) {
    std::cout << "\nApplying plane " << (i+1) << "/" << planes.size() << "..." << std::endl;

    auto result = clipPolyhedronByPlane(currentVertices, currentFaces,
                                       planes[i].first, planes[i].second);

    std::cout << "  Before: " << currentVertices.size() << " vertices, "
              << currentFaces.size() << " faces" << std::endl;
    std::cout << "  After:  " << result.vertices.size() << " vertices, "
              << result.faces.size() << " faces" << std::endl;

    if (result.fullyClipped) {
      std::cout << "  Polyhedron fully clipped!" << std::endl;
      break;
    }

    // Update for next iteration
    currentVertices = result.vertices;
    currentFaces = result.faces;
  }

  std::cout << "\nFinal result after all planes:" << std::endl;
  std::cout << "  Vertices: " << currentVertices.size() << std::endl;
  std::cout << "  Faces: " << currentFaces.size() << std::endl;

  // The final result should be the octant x>0, y>0, z>0
  // Original cube had one vertex at (1,1,1) which should remain
  bool hasPositiveOctantVertex = false;
  for (const auto& v : currentVertices) {
    if (v.x == 1 && v.y == 1 && v.z == 1) {
      hasPositiveOctantVertex = true;
      break;
    }
  }

  std::cout << "Contains (1,1,1) vertex: " << hasPositiveOctantVertex
            << " (expect 1)" << std::endl;
  assert(hasPositiveOctantVertex);

  std::cout << "PASS" << std::endl;
}

//------------------------------------------------------------------------------
// Main test driver
//------------------------------------------------------------------------------
int main() {
  std::cout << "Testing 3D Half-Space Clipping Layers" << std::endl;
  std::cout << "======================================" << std::endl;

  try {
    testPointClassification();
    testEdgeClipping();
    testFaceClipping();
    testPolyhedronClipping();
    testIterativeClipping();

    std::cout << "\n======================================" << std::endl;
    std::cout << "ALL TESTS PASSED" << std::endl;
    std::cout << "======================================" << std::endl;

    return 0;

  } catch (const std::exception& e) {
    std::cerr << "\nTEST FAILED: " << e.what() << std::endl;
    return 1;
  }
}
