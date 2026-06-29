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
  Point2<CoordType> gen0, gen1;
  Point2<CoordType> normalRay;

  // Output parameters
  shapes::BoxSide firstSide = shapes::BoxSide::L;
  shapes::BoxSide secondSide = shapes::BoxSide::L;
  shapes::BoxSide curSide;
  Point2<CoordType> p0, p1;

  // Input and output parameters
  bool inf0 = false, inf1 = false;

  // Returns true if edge should be skipped entirely
  bool doClipping(const Quantizer<2>& Q) {
    p1 = midPoint(gen1, gen0);
    if (inf0 && inf1) {
      normalRay = outwardRay(gen1, gen0);
      clipInfiniteRay(p1, -normalRay, Q.minBound, Q.maxBound, p0, curSide);
      firstSide = curSide;
      clipInfiniteRay(p0, normalRay, Q.minBound, Q.maxBound, p1, curSide);
      secondSide = curSide;
      return false;
    }
    POLY_ASSERT2(!inf0, "Cannot have only inf0");
    // Check if entire ray is external
    bool extRay = isRayExternal(rp0, normalRay, Q);
    if (extRay) {
      // If so, skip this ray entirely
      return true;
    }
    // Check if endpoint is outside the bounding box
    // If it is, leave it as the midpoint
    bool bounds1 = (!inf1) ? Q.inQBounds(rp1) : true;
    bool validp1 = (!inf1 && bounds1);
    if (validp1) {
      p1 = rp1.template type_cast<CoordType>();
    }
    bool bounds0 = Q.inQBounds(rp0);
    if (!bounds0 && !bounds1) {
      return true;
    }
    // If the origin is outside the bounds, clip it
    if (!bounds0) {
      clipInfiniteRay(p1, -normalRay, Q.minBound, Q.maxBound, p0, curSide);
      firstSide = curSide;
      inf0 = true;
    } else {
      p0 = rp0.template type_cast<CoordType>();
    }
    // Endpoint is outside bounding box
    if (!validp1) {
      clipInfiniteRay(p0, normalRay, Q.minBound, Q.maxBound, p1, curSide);
      secondSide = curSide;
      inf1 = true;
    }
    return false;
  }
};

}
#endif
