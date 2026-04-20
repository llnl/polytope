//------------------------------------------------------------------------
// QuantizedTessellation3d
//
// An intermediate representation for 3D tessellations in integer
// coordinates.
//------------------------------------------------------------------------
#ifndef __Polytope_QuantizedTessellation3d__
#define __Polytope_QuantizedTessellation3d__

#include <vector>
#include <utility>   // For std::pair

#include "Point.hh"
#include "polytope_tessellator_utilities.hh"

namespace polytope {

template<typename IntType, typename RealType>
struct QuantizedTessellation3d {
  typedef Point3<IntType> IntPoint;
  typedef Point3<RealType> RealPoint;
  
  RealType xmin[3], xmax[3], length, infRadius;
  std::vector<IntPoint> generators;

  std::vector<IntPoint> nodes;
  std::vector<std::pair<int, int> > edges;
  std::vector<std::vector<int> > faces;
  std::vector<std::vector<int> > cellFaces;

  // Construct with the given generators.  Finds the bounding limits and sets the
  // quantized generators.
  QuantizedTessellation3d(const std::vector<RealType>& points,
                          const std::vector<RealType>& boundaryPoints);

  // Construct with the given generators using the specified bounds.
  QuantizedTessellation3d(const std::vector<RealType>& points,
                          const RealType xmin_in[3],
                          const RealType xmax_in[3]);

  // Convert real coordinates to integers.
  void quantize(const RealType* realcoords, IntType* intcoords) const;

  // Convert int coordinates to reals.
  void dequantize(const IntType* intcoords, RealType* realcoords) const;

  // Read out the current QuantizedTessellation to regular Tessellation.
  void fillTessellation(Tessellation<3, RealType>& mesh) const;

  // Static properties about our coordinates.
  static IntType coordMin, coordMax;

private:
  // Internal construction method.
  void construct(const std::vector<RealType>& points);
};

}

#endif
