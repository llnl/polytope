//------------------------------------------------------------------------
// 2D clipping logic
//
//------------------------------------------------------------------------
#ifndef __Polytope_Clipping2D__
#define __Polytope_Clipping2D__

#include "EdgeUtils.hh"
#include "Intersections.hh"

namespace polytope {

// Raw geometry for a directed Voronoi primitive. The primitive is directed
// from start toward end. A missing endpoint denotes an infinite end.
template<typename CoordType>
struct DirectedVoronoiEdge {
public:
  bool hasStart = false;
  bool hasEnd = false;
  Point2<double> start;
  Point2<double> end;
  Point2<CoordType> bisectorPoint;
  Point2<CoordType> direction;

  void flipEdge() {
    if (hasStart) {
      end = start;
    }
    if (hasEnd) {
      start = end;
    }
    hasStart = !hasStart;
    hasEnd = !hasEnd;
    direction = -direction;
  }
};

// Clip directed Voronoi geometry to the quantization box. The assembler
// converts p0 and p1 to result.nodes indices after a successful clip.
template<typename CoordType>
class Clip2D {
public:
  explicit Clip2D(const DirectedVoronoiEdge<CoordType>& edge):
    mEdge(edge) {}

  // Output geometry. The assembler converts these to result.nodes indices.
  Point2<CoordType> p0;
  Point2<CoordType> p1;
  int startSide = -1;
  int endSide = -1;

  // Return true when the primitive intersects the quantization box.
  bool clip() {
    const auto& Q = Quantizer<2>::instance();
    startSide = -1;
    endSide = -1;

    const bool in0 = mEdge.hasStart && Q.inQBounds(mEdge.start);
    const bool in1 = mEdge.hasEnd && Q.inQBounds(mEdge.end);

    // No clipping is required for a finite edge wholly inside the box.
    if (in0 && in1) {
      p0 = round<2, CoordType>(mEdge.start);
      p1 = round<2, CoordType>(mEdge.end);
      return p0 != p1;
    }

    // Reject an out-of-box finite endpoint whose directed extension moves
    // farther away from the clipping box.
    if ((mEdge.hasStart && !in0 &&
         isRayExternal(mEdge.start, mEdge.direction)) ||
        (mEdge.hasEnd && !in1 &&
         isRayExternal(mEdge.end, -mEdge.direction))) {
      return false;
    }

    // The primitive is directed from p0 to p1. Replace a missing or
    // out-of-box endpoint with the corresponding clipping-box intersection.
    if (in0) {
      p0 = round<2, CoordType>(mEdge.start);
    } else {
      startSide =
        clipInfiniteRay(mEdge.bisectorPoint, -mEdge.direction, p0);
    }

    if (in1) {
      p1 = round<2, CoordType>(mEdge.end);
    } else {
      endSide =
        clipInfiniteRay(mEdge.bisectorPoint, mEdge.direction, p1);
    }

    return p0 != p1;
  }

private:
  const DirectedVoronoiEdge<CoordType>& mEdge;
};

}
#endif
