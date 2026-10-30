//------------------------------------------------------------------------
// 2D clipping logic
//
//------------------------------------------------------------------------
#ifndef __Polytope_Clipping2D__
#define __Polytope_Clipping2D__

#include "EdgeUtils.hh"
#include "Intersections.hh"

#include <iostream>

namespace polytope {

// Raw geometry for a directed Voronoi primitive. The primitive is directed
// from start toward end. A missing endpoint denotes an infinite end.
template<typename CoordType>
struct DirectedVoronoiEdge {
public:
  // Which generator is to the left of edge
  int leftGen = -1;
  bool hasStart = false;
  bool hasEnd = false;
  Point2<double> start;
  Point2<double> end;
  Point2<CoordType> bisectorPoint;
  Point2<CoordType> direction;

  void flipEdge() {
    if (hasStart && hasEnd) {
      Point2<double> temp = end;
      end = start;
      start = temp;
    } else if (hasStart) {
      end = start;
    } else if (hasEnd) {
      start = end;
    }
    hasStart = !hasStart;
    hasEnd = !hasEnd;
    direction = -direction;
  }
};

template<typename CoordType>
std::ostream&
operator<<(std::ostream& os,
           const DirectedVoronoiEdge<CoordType>& edge) {
  if (edge.hasStart) {
    os << "start " << edge.start << " ";
  }
  if (edge.hasEnd) {
    os << "end " << edge.end << " ";
  }
  if (edge.leftGen >= 0) {
    os << "leftgen " << edge.leftGen << " ";
  }
  os << "bisector " << edge.bisectorPoint << " dir " << edge.direction << std::endl;
  return os;
}

// Clip directed Voronoi geometry to the quantization box. The assembler
// converts p0 and p1 to result.nodes indices after a successful clip.
template<typename CoordType>
class Clip2D {
public:
  Clip2D(const DirectedVoronoiEdge<CoordType>& edge,
         const Quantizer<2>& Q):
    mEdge(edge), m_Q(Q) {}

  // Output geometry. The assembler converts these to result.nodes indices.
  Point2<CoordType> p0;
  Point2<CoordType> p1;
  int startSide = -1;
  int endSide = -1;

  // Return true when the primitive intersects the quantization box.
  bool clip() {
    startSide = -1;
    endSide = -1;

    // Check if direction is zero, implying degenerate point
    if (mEdge.direction.iszero()) {
      return false;
    }

    const bool in0 = mEdge.hasStart && m_Q.inQBounds(mEdge.start);
    const bool in1 = mEdge.hasEnd && m_Q.inQBounds(mEdge.end);

    // No clipping is required for a finite edge wholly inside the box.
    if (in0 && in1) {
      p0 = round<2, CoordType>(mEdge.start);
      p1 = round<2, CoordType>(mEdge.end);
      return p0 != p1;
    }

    // Reject an out-of-box finite endpoint whose directed extension moves
    // farther away from the clipping box.
    if ((mEdge.hasStart && !in0 &&
         isRayExternal(mEdge.start, mEdge.direction, m_Q)) ||
        (mEdge.hasEnd && !in1 &&
         isRayExternal(mEdge.end, -mEdge.direction, m_Q))) {
      return false;
    }

    // The primitive is directed from p0 to p1. Replace a missing or
    // out-of-box endpoint with the corresponding clipping-box intersection.
    if (in0) {
      p0 = round<2, CoordType>(mEdge.start);
    } else {
      startSide =
        clipInfiniteRay(mEdge.bisectorPoint, -mEdge.direction, p0, m_Q);
    }

    if (in1) {
      p1 = round<2, CoordType>(mEdge.end);
    } else {
      endSide =
        clipInfiniteRay(mEdge.bisectorPoint, mEdge.direction, p1, m_Q);
    }

    return p0 != p1;
  }

private:
  const DirectedVoronoiEdge<CoordType>& mEdge;
  const Quantizer<2>& m_Q;
};

}
#endif
