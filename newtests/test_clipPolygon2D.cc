// Comprehensive unit tests for 2D polygon clipping
//
// Tests the quantized 2D polygon clipping functionality including:
//   - Point-line classification (Layer 1)
//   - Edge-line clipping (Layer 2)
//   - Polygon-line clipping via Sutherland-Hodgman (Layer 3)
//   - Various polygon configurations and edge cases

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

#include "polytope.hh"
#include "Intersections.hh"
#include "Point.hh"
#include "polytope_test_utilities.hh"
#include "Boundary2D.hh"
#include "QuantTessellation.hh"
#include "BoostTessellator.hh"
#include "TriangleTessellator.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace polytope;
using namespace std;

namespace {

using CoordType = typename HashKey<2>::IntType;;
using IntPoint = Point2<CoordType>;

//------------------------------------------------------------------------------
// Test 1: Point-line classification (Layer 1)
//------------------------------------------------------------------------------
void testPointLineClassification(const int tnum) {
  cout << "\n=== Test " << tnum << ": Point-Line Classification ===" << endl;

  // Horizontal line from (0, 5) to (10, 5)
  IntPoint lineStart(0, 5);
  IntPoint lineEnd(10, 5);

  // Point above (left of line direction)
  IntPoint above(5, 10);
  int side1 = classifyPointByLine(above, lineStart, lineEnd);
  POLY_CHECK2(side1 == 1, "Point above line should return +1, got " << side1);

  // Point below (right of line direction)
  IntPoint below(5, 0);
  int side2 = classifyPointByLine(below, lineStart, lineEnd);
  POLY_CHECK2(side2 == -1, "Point below line should return -1, got " << side2);

  // Point on line
  IntPoint onLine(5, 5);
  int side3 = classifyPointByLine(onLine, lineStart, lineEnd);
  POLY_CHECK2(side3 == 0, "Point on line should return 0, got " << side3);

  // Vertical line from (5, 0) to (5, 10)
  IntPoint vLineStart(5, 0);
  IntPoint vLineEnd(5, 10);

  // Point to the left
  IntPoint left(0, 5);
  int side4 = classifyPointByLine(left, vLineStart, vLineEnd);
  POLY_CHECK2(side4 == 1, "Point left of vertical line should return +1, got " << side4);

  // Point to the right
  IntPoint right(10, 5);
  int side5 = classifyPointByLine(right, vLineStart, vLineEnd);
  POLY_CHECK2(side5 == -1, "Point right of vertical line should return -1, got " << side5);

  cout << "PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 2: Edge-line clipping (Layer 2)
//------------------------------------------------------------------------------
void testEdgeLineClipping(const int tnum) {
  cout << "\n=== Test " << tnum << ": Edge-Line Clipping ===" << endl;

  // Horizontal line from (0, 5) to (10, 5)
  IntPoint lineStart(0, 5);
  IntPoint lineEnd(10, 5);

  // Edge that crosses the line
  IntPoint e1Start(5, 0);
  IntPoint e1End(5, 10);
  IntPoint intersection1;
  int result1 = clipEdgeByLine(e1Start, e1End, lineStart, lineEnd, intersection1);
  POLY_CHECK2(result1 == 1, "Edge should intersect line, got result " << result1);
  POLY_CHECK2(intersection1.x == 5 && intersection1.y == 5,
              "Intersection should be at (5, 5), got (" << intersection1.x << ", " << intersection1.y << ")");

  // Edge that doesn't cross (parallel)
  IntPoint e2Start(0, 10);
  IntPoint e2End(10, 10);
  IntPoint intersection2;
  int result2 = clipEdgeByLine(e2Start, e2End, lineStart, lineEnd, intersection2);
  POLY_CHECK2(result2 == -1, "Parallel edge should return -1, got " << result2);

  // Edge with start on line
  IntPoint e3Start(5, 5);
  IntPoint e3End(5, 10);
  IntPoint intersection3;
  int result3 = clipEdgeByLine(e3Start, e3End, lineStart, lineEnd, intersection3);
  POLY_CHECK2(result3 == 0, "Edge starting on line should return 0, got " << result3);
  POLY_CHECK2(intersection3.x == 5 && intersection3.y == 5,
              "Intersection should be start point (5, 5)");

  cout << "PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 3: Clip square by horizontal line - keep bottom half
//------------------------------------------------------------------------------
void testClipSquareHorizontal(const int tnum) {
  cout << "\n=== Test " << tnum << ": Clip Square by Horizontal Line ===" << endl;

  // Unit square [0, 10] x [0, 10]
  vector<IntPoint> square = {
    IntPoint(0, 0),
    IntPoint(10, 0),
    IntPoint(10, 10),
    IntPoint(0, 10)
  };

  // Clip by horizontal line y = 5 (keep bottom half, which is left/above the line going right)
  IntPoint lineStart(0, 5);
  IntPoint lineEnd(10, 5);

  auto result = clipPolygonByLine(square, lineStart, lineEnd, true);

  POLY_CHECK2(!result.fullyClipped, "Square should not be fully clipped");
  POLY_CHECK2(!result.fullyRetained, "Square should be partially clipped");
  POLY_CHECK2(result.vertices.size() == 4,
              "Clipped square should have 4 vertices, got " << result.vertices.size());

  // Check expected vertices: (0, 0), (10, 0), (10, 5), (0, 5)
  bool hasBottomLeft = false, hasBottomRight = false;
  bool hasTopRight = false, hasTopLeft = false;

  for (const auto& v : result.vertices) {
    if (v.x == 0 && v.y == 0) hasBottomLeft = true;
    if (v.x == 10 && v.y == 0) hasBottomRight = true;
    if (v.x == 10 && v.y == 5) hasTopRight = true;
    if (v.x == 0 && v.y == 5) hasTopLeft = true;
  }

  POLY_CHECK2(hasBottomLeft && hasBottomRight && hasTopRight && hasTopLeft,
              "Clipped square should have correct vertices");

  cout << "PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 4: Clip square by vertical line - keep left half
//------------------------------------------------------------------------------
void testClipSquareVertical(const int tnum) {
  cout << "\n=== Test " << tnum << ": Clip Square by Vertical Line ===" << endl;

  // Unit square [0, 10] x [0, 10]
  vector<IntPoint> square = {
    IntPoint(0, 0),
    IntPoint(10, 0),
    IntPoint(10, 10),
    IntPoint(0, 10)
  };

  // Clip by vertical line x = 5 (going upward, keep left side)
  IntPoint lineStart(5, 0);
  IntPoint lineEnd(5, 10);

  auto result = clipPolygonByLine(square, lineStart, lineEnd, true);

  POLY_CHECK2(!result.fullyClipped, "Square should not be fully clipped");
  POLY_CHECK2(!result.fullyRetained, "Square should be partially clipped");
  POLY_CHECK2(result.vertices.size() == 4,
              "Clipped square should have 4 vertices, got " << result.vertices.size());

  cout << "PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 5: Clip square by diagonal line
//------------------------------------------------------------------------------
void testClipSquareDiagonal(const int tnum) {
  cout << "\n=== Test " << tnum << ": Clip Square by Diagonal Line ===" << endl;

  // Unit square [0, 10] x [0, 10]
  vector<IntPoint> square = {
    IntPoint(0, 0),
    IntPoint(10, 0),
    IntPoint(10, 10),
    IntPoint(0, 10)
  };

  // Clip by diagonal line from (0, 0) to (10, 10) - keep lower-left side
  IntPoint lineStart(0, 0);
  IntPoint lineEnd(10, 10);

  auto result = clipPolygonByLine(square, lineStart, lineEnd, true);

  POLY_CHECK2(!result.fullyClipped, "Square should not be fully clipped");
  POLY_CHECK2(!result.fullyRetained, "Square should be partially clipped");
  POLY_CHECK2(result.vertices.size() == 3,
              "Clipped square should form triangle with 3 vertices, got " << result.vertices.size());

  cout << "PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 6: Clip triangle - partial intersection
//------------------------------------------------------------------------------
void testClipTriangle(const int tnum) {
  cout << "\n=== Test " << tnum << ": Clip Triangle ===" << endl;

  // Triangle
  vector<IntPoint> triangle = {
    IntPoint(5, 0),
    IntPoint(10, 10),
    IntPoint(0, 10)
  };

  // Clip by horizontal line y = 5
  IntPoint lineStart(0, 5);
  IntPoint lineEnd(10, 5);

  auto result = clipPolygonByLine(triangle, lineStart, lineEnd, true);

  POLY_CHECK2(!result.fullyClipped, "Triangle should not be fully clipped");
  POLY_CHECK2(!result.fullyRetained, "Triangle should be partially clipped");
  POLY_CHECK2(result.vertices.size() == 4,
              "Clipped triangle should have 4 vertices, got " << result.vertices.size());

  cout << "PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 7: Polygon fully inside (fully retained)
//------------------------------------------------------------------------------
void testPolygonFullyInside(const int tnum) {
  cout << "\n=== Test " << tnum << ": Polygon Fully Inside ===" << endl;

  // Small square [2, 8] x [2, 8]
  vector<IntPoint> square = {
    IntPoint(2, 2),
    IntPoint(8, 2),
    IntPoint(8, 8),
    IntPoint(2, 8)
  };

  // Clip by line far away: y = 0 (horizontal line below square)
  IntPoint lineStart(0, 0);
  IntPoint lineEnd(10, 0);

  auto result = clipPolygonByLine(square, lineStart, lineEnd, true);

  POLY_CHECK2(!result.fullyClipped, "Square should not be clipped");
  POLY_CHECK2(result.fullyRetained, "Square should be fully retained");
  POLY_CHECK2(result.vertices.size() == 4,
              "Square should keep all 4 vertices, got " << result.vertices.size());

  cout << "PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 8: Polygon fully outside (fully clipped)
//------------------------------------------------------------------------------
void testPolygonFullyOutside(const int tnum) {
  cout << "\n=== Test " << tnum << ": Polygon Fully Outside ===" << endl;

  // Square [0, 10] x [0, 10]
  vector<IntPoint> square = {
    IntPoint(0, 0),
    IntPoint(10, 0),
    IntPoint(10, 10),
    IntPoint(0, 10)
  };

  // Clip by line far above: y = 20, keep below (which is left side going right)
  IntPoint lineStart(0, 20);
  IntPoint lineEnd(10, 20);

  auto result = clipPolygonByLine(square, lineStart, lineEnd, true);

  POLY_CHECK2(result.fullyClipped, "Square should be fully clipped");
  POLY_CHECK2(!result.fullyRetained, "Square should not be retained");
  POLY_CHECK2(result.vertices.size() == 0,
              "Fully clipped polygon should have 0 vertices, got " << result.vertices.size());

  cout << "PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 9: Polygon edge exactly on clipping line
//------------------------------------------------------------------------------
void testPolygonEdgeOnLine(const int tnum) {
  cout << "\n=== Test " << tnum << ": Polygon Edge on Clipping Line ===" << endl;

  // Square [0, 10] x [0, 10]
  vector<IntPoint> square = {
    IntPoint(0, 0),
    IntPoint(10, 0),
    IntPoint(10, 10),
    IntPoint(0, 10)
  };

  // Clip by line along bottom edge: y = 0
  IntPoint lineStart(0, 0);
  IntPoint lineEnd(10, 0);

  auto result = clipPolygonByLine(square, lineStart, lineEnd, true);

  POLY_CHECK2(!result.fullyClipped, "Square should not be fully clipped");
  // Should keep the square since vertices on the line are kept
  POLY_CHECK2(result.vertices.size() >= 3,
              "Square should have at least 3 vertices after clipping");

  cout << "PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 10: Sequential clipping (multiple edges)
//------------------------------------------------------------------------------
void testSequentialClipping(const int tnum) {
  cout << "\n=== Test " << tnum << ": Sequential Clipping ===" << endl;

  // Square [0, 10] x [0, 10]
  vector<IntPoint> square = {
    IntPoint(0, 0),
    IntPoint(10, 0),
    IntPoint(10, 10),
    IntPoint(0, 10)
  };

  // Clip by bottom edge: y = 2, keep above (left side)
  IntPoint line1Start(0, 2);
  IntPoint line1End(10, 2);
  auto result1 = clipPolygonByLine(square, line1Start, line1End, true);

  POLY_CHECK2(!result1.fullyClipped, "First clip should not fully clip");

  // Clip by left edge: x = 2, keep right (left side of line going up)
  IntPoint line2Start(2, 0);
  IntPoint line2End(2, 10);
  auto result2 = clipPolygonByLine(result1.vertices, line2Start, line2End, true);

  POLY_CHECK2(!result2.fullyClipped, "Second clip should not fully clip");

  // Clip by top edge: y = 8, keep below (left side of line going right)
  IntPoint line3Start(0, 8);
  IntPoint line3End(10, 8);
  auto result3 = clipPolygonByLine(result2.vertices, line3Start, line3End, true);

  POLY_CHECK2(!result3.fullyClipped, "Third clip should not fully clip");

  // Clip by right edge: x = 8, keep left (left side of line going up)
  IntPoint line4Start(8, 0);
  IntPoint line4End(8, 10);
  auto result4 = clipPolygonByLine(result3.vertices, line4Start, line4End, true);

  POLY_CHECK2(!result4.fullyClipped, "Fourth clip should not fully clip");
  POLY_CHECK2(result4.vertices.size() == 4,
              "Final clipped square should have 4 vertices, got " << result4.vertices.size());

  // Should end up with [2, 8] x [2, 8] square
  // Check all vertices are within bounds
  for (const auto& v : result4.vertices) {
    POLY_CHECK2(v.x >= 2 && v.x <= 8 && v.y >= 2 && v.y <= 8,
                "Vertex should be in range [2, 8]");
  }

  cout << "PASS" << endl;
}

//------------------------------------------------------------------------------
// Test 10: Complex clipping (multiple edges)
//------------------------------------------------------------------------------
void testSquare(const int tnum, bool boostTess) {
  cout << "\n=== Test " << tnum << ": Square Clipping ===" << endl;
  std::string outname = (boostTess) ? "boost" : "triangle";
  outname += std::to_string(tnum);
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(0);
  Quantizer<2> Q(boundary.mQ);
  double lov = -0.5;
  double hiv = 0.5;
  double len = hiv - lov;
  unsigned Nx = 2;
  double dx = len/double(Nx);
  vector<double> points = {0.05, 0.025, 0.025, 0.05, 0.05, -0.05, -0.05, -0.05, -0.05, 0.05};
  QuantTessellation<2> quantMesh(Q, points);
  QuantPLC<2> QPLC(boundary.mPLC, Q, boundary.mPLCpoints);
  if (boostTess) {
    BoostTessellator boost(Q);
    boost.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, boost);
  } else {
    TriangleTessellator tri(Q);
    tri.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, tri);
  }
  Tessellation<2, double> mesh;
  quantMesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  outputMesh(mesh, outname, points, 0, 0.0);
}

//------------------------------------------------------------------------------
// Test 10: Complex clipping (multiple edges)
//------------------------------------------------------------------------------
void testSquareTriangle(const int tnum, bool boostTess) {
  cout << "\n=== Test " << tnum << ": Square Triangle Clipping ===" << endl;
  std::string outname = (boostTess) ? "boost" : "triangle";
  outname += std::to_string(tnum);
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(0);
  vector<double> points = {0.05, 0.025, 0.025, 0.05, 0.05, -0.05, -0.05, -0.05, -0.05, 0.05};
  // Add a triangle hole
  vector<double> newPoints = {0.3, -0.4, 0.2, 0.4, 0.2, -0.4};
  auto Nf = newPoints.size()/2;
  auto N = boundary.mPLCpoints.size()/2;
  copy(newPoints.begin(), newPoints.end(), back_inserter(boundary.mPLCpoints));
  boundary.mPLC.holes = vector<vector<vector<int>>>(1);
  boundary.mPLC.holes[0].resize(Nf);
  for (int i = 0; i < Nf; ++i) {
    unsigned fbegin = N + i;
    unsigned fend = N + (i+1)%Nf;
    boundary.mPLC.holes[0][i].resize(2);
    boundary.mPLC.holes[0][i][0] = fbegin;
    boundary.mPLC.holes[0][i][1] = fend;
  }
  //boundary.finalize();
  Quantizer<2> Q(boundary.mQ);
  QuantTessellation<2> quantMesh(Q, points);
  QuantPLC<2> QPLC(boundary.mPLC, Q, boundary.mPLCpoints);
  if (boostTess) {
    BoostTessellator boost(Q);
    boost.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, boost);
  } else {
    TriangleTessellator tri(Q);
    tri.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, tri);
  }
  Tessellation<2, double> mesh;
  quantMesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  outputMesh(mesh, outname, points, 0, 0.0);
}

//------------------------------------------------------------------------------
// Test 10: Complex clipping (multiple edges)
//------------------------------------------------------------------------------
void testDiamond(const int tnum, bool boostTess) {
  cout << "\n=== Test " << tnum << ": Complex Clipping ===" << endl;
  std::string outname = (boostTess) ? "boost" : "triangle";
  outname += std::to_string(tnum);
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(0);
  Quantizer<2> Q(boundary.mQ);
  double lov = -0.5;
  double hiv = 0.5;
  vector<double> points = {lov, 0, hiv, 0, 0, lov, 0, hiv};
  QuantTessellation<2> quantMesh(Q, points);
  QuantPLC<2> QPLC(boundary.mPLC, Q, boundary.mPLCpoints);
  if (boostTess) {
    BoostTessellator boost(Q);
    boost.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, boost);
  } else {
    TriangleTessellator tri(Q);
    tri.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, tri);
  }
  Tessellation<2, double> mesh;
  quantMesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  outputMesh(mesh, outname, points, 0, 0.0);
}


} // anonymous namespace

//------------------------------------------------------------------------------
// Main
//------------------------------------------------------------------------------
int main(int argc, char** argv) {
#ifdef POLYTOPE_ENABLE_MPI
  MPI_Init(&argc, &argv);
#endif

  try {
    int test = 1;
    bool boost = true;
    for (int i = 0; i < 2; ++i) {
      testSquare(test++, boost);
      testSquareTriangle(test++, boost);
      testDiamond(test++, boost);
      boost = false;
    }
    testPointLineClassification(test++);
    testEdgeLineClipping(test++);
    testClipSquareHorizontal(test++);
    testClipSquareVertical(test++);
    testClipSquareDiagonal(test++);
    testClipTriangle(test++);
    testPolygonFullyInside(test++);
    testPolygonFullyOutside(test++);
    testPolygonEdgeOnLine(test++);
    testSequentialClipping(test++);

    cout << "\n=== ALL TESTS PASSED ===" << endl;
  } catch (const exception& e) {
    cout << "\n=== TEST FAILED WITH EXCEPTION ===" << endl;
    cout << e.what() << endl;
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
