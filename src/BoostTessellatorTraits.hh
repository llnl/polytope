//------------------------------------------------------------------------
// BoostTessellatorTraits
//
// Traits class to help with integer and floating point operations in
// BoostTessellator. We have some freedom in how the tessellator is
// operated. Boost.Voronoi takes integer input (default is 32-bit) and
// returns floating-point Voronoi vertices. Unlike TriangleTessellator,
// in which quantizing overlapping triangle circumcenters is necessary for 
// resolving degenerate mesh configurations, Boost constructs a valid
// topology that resolves degeneracies on its own. Doing all operations
// in integers is not AS necessary for Boost as it is for Triangle.
//
// The following trait class allows you to specify a type choice for
// all internal operations in BoostTessellator. The trait class provides
// a simple interface for moving back and forth from floating point and 
// integers, regardless of whether you chose floating point data as the
// internal data type or integers.
//
// deboost: Take a floating-point vertex from the Boost.Voronoi struct
//          and return a Polytope Point of the location
//
// quantize: Hash a floating point vertex to the integer mesh.
//           If input is integer, just return what you handed in.
//
// dequantize: Return an integer (hashed) point to floating point value.
//             If input is floating point, just return what was handed in.
//
// project: Project a point to the infinite bounding sphere using a vertex
//          location (int or float point) and a floating point direction.
//------------------------------------------------------------------------
#ifndef __Polytope_BoostTessellatorTraits__
#define __Polytope_BoostTessellatorTraits__
#include "Point.hh"
#include "QuantizedCoordinates.hh"

// The abstraction
template<typename RealType, typename BaseType> struct BoostTessellatorTraits {};

//------------------------------------------------------------------------
// Floating Point Specialization
//
// Base type is floating point. The data out of Boost has to be rescaled based on
// the original quantization of the generators. The return type is a RealPoint.
// Calling quantize should do nothing. We're not working in the ints with this
// option turned on.
//------------------------------------------------------------------------
template<typename RealType>
struct BoostTessellatorTraits<RealType, RealType> {
  typedef typename polytope::DimensionTraits<2, RealType>::CoordHash IntType;
  typedef RealType CoordType;
  typedef polytope::Point2<RealType> PointType;
  typedef polytope::Point2<RealType> RealPoint;
  typedef polytope::Point2<IntType>  IntPoint;
  static PointType deboost(const polytope::QuantizedCoordinates<2, RealType>& coords, 
                           const RealType* pointIn) {
    return RealPoint(coords.low_inner.x + coords.delta*pointIn[0],
                     coords.low_inner.y + coords.delta*pointIn[1]);
  }
  static PointType quantize(const polytope::QuantizedCoordinates<2, RealType>& coords, 
                            const RealType* pointIn) {
    return PointType(pointIn[0], pointIn[1]);
  }
  static RealPoint dequantize(const polytope::QuantizedCoordinates<2, RealType>& coords, 
                              const PointType pointIn) {
    return pointIn;
  }
  static PointType project(const polytope::QuantizedCoordinates<2, RealType>& coords, 
                           const RealType* endPoint,
                           const RealType* direction) {
    return coords.projectPoint(endPoint, direction);
  }
};


//------------------------------------------------------------------------
// Integer Specialization
//
// Base type is 64-bit integer. The data out of Boost has to be floored. Although
// it's scaled to be on the integer grid, it's actually floating point valued.
// The return type is an IntPoint. Calling quantize should snap a floating point
// value to the integer grid. Geometric operation use the ints with this option on.
//------------------------------------------------------------------------
template<typename RealType>
struct BoostTessellatorTraits<RealType, typename polytope::DimensionTraits<2, RealType>::CoordHash> {
  typedef typename polytope::DimensionTraits<2, RealType>::CoordHash IntType;
  typedef IntType CoordType;
  typedef polytope::Point2<IntType>  PointType;
  typedef polytope::Point2<RealType> RealPoint;
  typedef polytope::Point2<IntType>  IntPoint;
  static PointType deboost(const polytope::QuantizedCoordinates<2, RealType>& coords, 
                           const RealType* pointIn) {
    return PointType(floor(pointIn[0]), floor(pointIn[1]));
  }
  static PointType quantize(const polytope::QuantizedCoordinates<2, RealType>& coords, 
                            const RealType* pointIn) {
    return coords.quantize(pointIn);
  }
  static RealPoint dequantize(const polytope::QuantizedCoordinates<2, RealType>& coords, 
                              const PointType pointIn) {
    return coords.dequantize(&pointIn.x);
  }
  static PointType project(const polytope::QuantizedCoordinates<2, RealType>& coords, 
                           const IntType* endPoint,
                           const RealType* direction) {
    const IntPoint iep   = IntPoint(endPoint[0], endPoint[1]);
    const RealPoint tmp  = coords.dequantize(&iep.x);
    const RealPoint pinf = coords.projectPoint(&tmp.x, direction);
    return coords.quantize(&pinf.x);
  }
};
#endif
