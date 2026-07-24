//------------------------------------------------------------------------
// 2D clipping logic
//
//------------------------------------------------------------------------
#ifndef __Polytope_Clipping2D__
#define __Polytope_Clipping2D__

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
  Point2<CoordType> p0, p1;

  // Input and output parameters
  bool inf0 = false, inf1 = false;

  // Returns true if edge should be skipped entirely
  bool doClipping() {
    auto& Q = Quantizer<2>::instance();
    bool clip1 = true;
    bool clip2 = true;
    Point2<CoordType> m = midPoint(gen1, gen0);
    if (inf0 && inf1) {
      clip1 = clipInfiniteRay(m, -normalRay, p0, firstSide);
      clip2 = clipInfiniteRay(m, normalRay, p1, secondSide);
      POLY_ASSERT(clip1 == clip2);
      return !clip1;
    }
    // if (!inf0 && !inf1) {
    //   Point2<double> a, b;
    //   bool clippedFirst = false;
    //   bool clippedSecond = false;
    //   if (!clipFiniteSegmentToBox<CoordType>(rp0, rp1, a, b, firstSide, secondSide,
    //                                          clippedFirst, clippedSecond)) {
    //     return true;
    //   }
    //   p0 = round<2, CoordType>(a+0.5);
    //   p1 = round<2, CoordType>(b+0.5);
    //   inf0 = clippedFirst;
    //   inf1 = clippedSecond;
    //   return false;
    // }
    bool bounds0 = (!inf0) ? Q.inQBounds(rp0) : true;
    bool validp0 = (!inf0 && bounds0);
    if (validp0) {
      p0 = rp0.template type_cast<CoordType>();
    }
    bool bounds1 = (!inf1) ? Q.inQBounds(rp1) : true;
    bool validp1 = (!inf1 && bounds1);
    if (validp1) {
      p1 = rp1.template type_cast<CoordType>();
    }
    if (validp0 && validp1) {
      return false;
    }
    if ((!bounds0 && isRayExternal(rp0, normalRay)) ||
        (!bounds1 && isRayExternal(rp1, -normalRay))) {
      return true;
    }
    if (!validp0) {
      clip1 = clipInfiniteRay(m, -normalRay, p0, firstSide);
      inf0 = true;
      if (!clip1) return true;
    }
    if (!validp1) {
      clip2 = clipInfiniteRay(m, normalRay, p1, secondSide);
      inf1 = true;
      if (!clip2) return true;
    }
    return false;
  }
};

}
#endif
