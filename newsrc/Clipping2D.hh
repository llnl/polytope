//------------------------------------------------------------------------
// 2D clipping logic
//
//------------------------------------------------------------------------
#ifndef POLYTOPE_CLIPPING2D_HH
#define POLYTOPE_CLIPPING2D_HH

#include "Intersections.hh"

namespace polytope {

template<typename CoordType>
struct Clip2D {
  // Input parameters
  Point2<double> rp0, rp1;
  // Generator nodes
  Point2<CoordType> gen0, gen1;

  // Output parameters
  shapes::BoxSide firstSide; // In case both are infinite
  shapes::BoxSide curSide;
  Point2<CoordType> p0, p1;
  bool bothInf = false;
  bool bounds0 = true, bounds1 = true;

  // Input and output parameters
  bool inf0 = false, inf1 = false;

  // Returns true if edge should be skipped entirely
  bool doClipping(const Quantizer<2>& Q) {
    if (!inf0) {
      bounds0 = Q.inQBounds(rp0);
    }
    if (!inf1) {
      bounds1 = Q.inQBounds(rp1);
    }
    // If points are near boundary edges, skip point
    if ((!bounds0 && inf1) || (!bounds1 && inf0) || (!bounds0 && !bounds1)) {
      return true;
    }
    bool doReturn = false;
    if (bounds0 && !inf0) {
      p0 = round<2, CoordType>(rp0);
      doReturn = true;
    }
    if (bounds1 && !inf1) {
      p1 = round<2, CoordType>(rp1);
      if (doReturn) return false;
    }
    // If both vertices are infinite, start ray at the midpoint between the generators
    if (inf0 && inf1) {
      p1 = midPoint(gen1, gen0);
      bothInf = true;
    }
    if (inf0 || !bounds0) {
      Point2<CoordType> outwardRay = normalRay(gen1, gen0);
      clipInfiniteRay(p1, outwardRay, Q.minBound, Q.maxBound, p0, curSide);
      if (bothInf) {
        firstSide = curSide;
      }
      inf0 = true;
    }
    if (inf1 || !bounds1) {
      Point2<CoordType> outwardRay = normalRay(gen0, gen1);
      clipInfiniteRay(p0, outwardRay, Q.minBound, Q.maxBound, p1, curSide);
      inf1 = true;
    }
    return false;
  }
};

}
#endif
