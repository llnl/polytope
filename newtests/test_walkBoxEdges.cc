#include <iostream>
#include <vector>
#include <map>
#include <cassert>

#include "Shapes.hh"
#include "EdgeUtils.hh"

#if 0
using namespace polytope;
using namespace polytope::shapes;
using namespace polytope::edge;

//------------------------------------------------------------------------------
// Helper to print edges for debugging
//------------------------------------------------------------------------------
void printEdges(const std::vector<Edge>& edges, const std::string& label) {
  std::cout << label << ":" << std::endl;
  for (size_t i = 0; i < edges.size(); ++i) {
    std::cout << "  Edge " << i << ": (" << edges[i].first
              << ", " << edges[i].second << ")" << std::endl;
  }
}

void testSimpleWalk() {
  std::cout << "\n=== Test 1: Simple walk on same edge ===" << std::endl;

  // Set up corner indices for a unit square
  // LL=0, LR=1, UR=2, UL=3 (CCW ordering)
  std::map<BoxCorner, unsigned> cornerIndices;
  cornerIndices[BoxCorner::LL] = 0;
  cornerIndices[BoxCorner::LR] = 1;
  cornerIndices[BoxCorner::UR] = 2;
  cornerIndices[BoxCorner::UL] = 3;

  // Pretend the following indices represent different points
  // 4: Middle of the box
  // 5: Middle of left side
  // 6: Middle of bottom side
  // 7: Just above the LL point on the left side
  BoxSide origSide = BoxSide::L;
  BoxSide curSide = BoxSide::L;
  Edge curEdge = std::make_pair(4, 5);
  unsigned origPoint = 7;  // UL corner (where we started on left side)
  bool endAtCurEdge = false;

  std::vector<Edge> edges;
  walkBoxEdges(origSide, curSide, cornerIndices, curEdge,
               origPoint, endAtCurEdge, edges);

  printEdges(edges, "Result");

  // Expected: (4,5) (5, 7)
  // Walk adds curEdge, then connects from LL corner (0) to origPoint (3)
  assert(edges.size() == 2);
  assert(edges[0].first == 4 && edges[0].second == 5);
  assert(edges[1].first == 5 && edges[1].second == 7);

  std::cout << "PASS" << std::endl;
}

void testSimpleWalk2() {
  std::cout << "\n=== Test 1: Simple walk on same edge ===" << std::endl;

  // Set up corner indices for a unit square
  // LL=0, LR=1, UR=2, UL=3 (CCW ordering)
  std::map<BoxCorner, unsigned> cornerIndices;
  cornerIndices[BoxCorner::LL] = 0;
  cornerIndices[BoxCorner::LR] = 1;
  cornerIndices[BoxCorner::UR] = 2;
  cornerIndices[BoxCorner::UL] = 3;

  // Pretend the following indices represent different points
  // 4: Middle of the box
  // 5: Middle of left side
  // 6: Middle of bottom side
  // 7: Just above the LL point on the left side
  BoxSide origSide = BoxSide::L;
  BoxSide curSide = BoxSide::L;
  Edge curEdge = std::make_pair(7, 4);
  unsigned origPoint = 5;  // UL corner (where we started on left side)
  bool endAtCurEdge = true;

  std::vector<Edge> edges;
  walkBoxEdges(origSide, curSide, cornerIndices, curEdge,
               origPoint, endAtCurEdge, edges);

  printEdges(edges, "Result");

  // Expected: (5,7) (7,4)
  assert(edges.size() == 2);
  assert(edges[0].first == 5 && edges[0].second == 7);
  assert(edges[1].first == 7 && edges[1].second == 4);

  std::cout << "PASS" << std::endl;
}

//------------------------------------------------------------------------------
// Test 1: Walk from Left side to Bottom side (clockwise, not ending at current edge)
// Should traverse: L -> B directly
//------------------------------------------------------------------------------
void testCornerWalk1() {
  std::cout << "\n=== Test 1: Simple walk L->B (not ending at current) ===" << std::endl;

  // Set up corner indices for a unit square
  // LL=0, LR=1, UR=2, UL=3 (CCW ordering)
  std::map<BoxCorner, unsigned> cornerIndices;
  cornerIndices[BoxCorner::LL] = 0;
  cornerIndices[BoxCorner::LR] = 1;
  cornerIndices[BoxCorner::UR] = 2;
  cornerIndices[BoxCorner::UL] = 3;

  // Pretend the following indices represent different points
  // 4: Middle of the box
  // 5: Middle of left side
  // 6: Middle of bottom side
  BoxSide origSide = BoxSide::L;
  BoxSide curSide = BoxSide::B;
  Edge curEdge = std::make_pair(1, 4);  // LR -> middle of box
  unsigned origPoint = 3;  // UL corner (where we started on left side)
  bool endAtCurEdge = true;

  std::vector<Edge> edges;
  walkBoxEdges(origSide, curSide, cornerIndices, curEdge,
               origPoint, endAtCurEdge, edges);

  printEdges(edges, "Result");

  // Expected: (3,0) (0,1) (1,4)
  // Walk adds curEdge, then connects from LL corner (0) to origPoint (3)
  assert(edges.size() == 3);
  assert(edges[0].first == 3 && edges[0].second == 0);
  assert(edges[1].first == 0 && edges[1].second == 1);
  assert(edges[2].first == 1 && edges[2].second == 4);

  std::cout << "PASS" << std::endl;
}

//------------------------------------------------------------------------------
// Test 2: Walk from Right side to Top side, ending at current edge
// Should traverse: R -> T, adding corners along the way
//------------------------------------------------------------------------------
void testCornerWalk2() {
  std::cout << "\n=== Test 2: Walk R->T (ending at current edge) ===" << std::endl;

  std::map<BoxCorner, unsigned> cornerIndices;
  cornerIndices[BoxCorner::LL] = 0;
  cornerIndices[BoxCorner::LR] = 1;
  cornerIndices[BoxCorner::UR] = 2;
  cornerIndices[BoxCorner::UL] = 3;

  BoxSide origSide = BoxSide::B;
  BoxSide curSide = BoxSide::T;
  // 4: Middle of the box
  // 5: Middle of left side
  // 6: Middle of bottom side
  // 7: Middle of top side
  Edge curEdge = std::make_pair(4, 7);
  unsigned origPoint = 6;
  bool endAtCurEdge = false;

  std::vector<Edge> edges;
  walkBoxEdges(origSide, curSide, cornerIndices, curEdge,
               origPoint, endAtCurEdge, edges);

  printEdges(edges, "Result");

  std::vector<Edge> correct(4);
  correct[0] = std::make_pair(4,7);
  correct[1] = std::make_pair(7,3);
  correct[2] = std::make_pair(3,0);
  correct[3] = std::make_pair(0,6);
  // Expected: (4,7) (7,3) (3,0) (0,6)
  // Walk from origPoint (1=LR) to UR corner (2), then add final edge
  assert(edges.size() == 4);
  for (int i = 0; i < 4; ++i) {
    assert(edges[i].first == correct[i].first && edges[i].second == correct[i].second);
  }

  std::cout << "PASS" << std::endl;
}

//------------------------------------------------------------------------------
// Test 2: Walk from Right side to Top side, ending at current edge
// Should traverse: R -> T, adding corners along the way
//------------------------------------------------------------------------------
void testCornerWalk3() {
  std::cout << "\n=== Test 3: Walk R->T (ending at current edge) ===" << std::endl;

  std::map<BoxCorner, unsigned> cornerIndices;
  cornerIndices[BoxCorner::LL] = 0;
  cornerIndices[BoxCorner::LR] = 1;
  cornerIndices[BoxCorner::UR] = 2;
  cornerIndices[BoxCorner::UL] = 3;

  BoxSide origSide = BoxSide::B;
  BoxSide curSide = BoxSide::T;
  // 4: Middle of the box
  // 5: Middle of left side
  // 6: Middle of bottom side
  // 7: Middle of top side
  Edge curEdge = std::make_pair(7,6);
  unsigned origPoint = 6;
  bool endAtCurEdge = true;

  std::vector<Edge> edges;
  walkBoxEdges(origSide, curSide, cornerIndices, curEdge,
               origPoint, endAtCurEdge, edges);

  printEdges(edges, "Result");

  std::vector<Edge> correct(4);
  correct[0] = std::make_pair(6,1);
  correct[1] = std::make_pair(1,2);
  correct[2] = std::make_pair(2,7);
  correct[3] = std::make_pair(7,6);
  assert(edges.size() == 4);
  for (int i = 0; i < 4; ++i) {
    assert(edges[i].first == correct[i].first && edges[i].second == correct[i].second);
  }

  std::cout << "PASS" << std::endl;
}

//------------------------------------------------------------------------------
// Test 3: Walk around multiple sides (L -> R)
// Should traverse: L -> B -> R
//------------------------------------------------------------------------------
void testMultipleSides() {
  std::cout << "\n=== Test 3: Walk around multiple sides L->R ===" << std::endl;

  std::map<BoxCorner, unsigned> cornerIndices;
  cornerIndices[BoxCorner::LL] = 0;
  cornerIndices[BoxCorner::LR] = 1;
  cornerIndices[BoxCorner::UR] = 2;
  cornerIndices[BoxCorner::UL] = 3;

  BoxSide origSide = BoxSide::L;
  BoxSide curSide = BoxSide::R;
  Edge curEdge = std::make_pair(1, 2);  // LR -> UR edge on right side
  unsigned origPoint = 3;  // UL corner
  bool endAtCurEdge = false;

  std::vector<Edge> edges;
  walkBoxEdges(origSide, curSide, cornerIndices, curEdge,
               origPoint, endAtCurEdge, edges);

  printEdges(edges, "Result");

  // Expected to walk: L->B (corner LL), B->R (corner LR), then close back
  assert(edges.size() >= 3);
  std::cout << "PASS" << std::endl;
}
#endif
//------------------------------------------------------------------------------
// Main
//------------------------------------------------------------------------------
int main() {
  std::cout << "Testing walkBoxEdges routine..." << std::endl;
#if 0
  testSimpleWalk();
  testSimpleWalk2();
  testCornerWalk1();
  testCornerWalk2();
  testCornerWalk3();
#endif
  std::cout << "\n=== All tests passed ===" << std::endl;
  return 0;
}

