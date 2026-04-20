//------------------------------------------------------------------------
// QuantizedTessellation3d
//
// An intermediate representation for 3D tessellations in integer
// coordinates.
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <map>

#include "QuantizedTessellation3d.hh"
#include "polytope.hh"
#include "polytope_tessellator_utilities.hh"

namespace polytope {
  
//------------------------------------------------------------------------------
// Construct with the given generators.  Finds the bounding limits and sets the
// quantized generators.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
QuantizedTessellation3d<IntType, RealType>::
QuantizedTessellation3d(const std::vector<RealType>& points,
                        const std::vector<RealType>& boundaryPoints) {
  geometry::computeBoundingBox<3, RealType>(&boundaryPoints[0], boundaryPoints.size(), true, xmin, xmax);
  length = std::max(xmax[0] - xmin[0], std::max(xmax[1] - xmin[1], xmax[2] - xmin[2]));
  xmin[0] -= 2.0*length;
  xmin[1] -= 2.0*length;
  xmin[2] -= 2.0*length;
  xmax[0] += 2.0*length;
  xmax[1] += 2.0*length;
  xmax[2] += 2.0*length;
  length *= 5.0;
  this->construct(points);
}

//------------------------------------------------------------------------------
// Construct with the given generators using the specified bounds.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
QuantizedTessellation3d<IntType, RealType>::
QuantizedTessellation3d(const std::vector<RealType>& points,
                        const RealType xmin_in[3],
                        const RealType xmax_in[3]) {
  xmin[0] = xmin_in[0];
  xmin[1] = xmin_in[1];
  xmin[2] = xmin_in[2];
  xmax[0] = xmax_in[0];
  xmax[1] = xmax_in[1];
  xmax[2] = xmax_in[2];
  length = std::max(xmax[0] - xmin[0], std::max(xmax[1] - xmin[1], xmax[2] - xmin[2]));
  this->construct(points);
}

//------------------------------------------------------------------------------
// Convert real coordinates to integers.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
void
QuantizedTessellation3d<IntType, RealType>::
quantize(const RealType* realcoords, IntType* intcoords) const {
  const RealType dx = length/(coordMax - coordMin);
  intcoords[0] = coordMin + IntType((realcoords[0] - xmin[0])/dx);
  intcoords[1] = coordMin + IntType((realcoords[1] - xmin[1])/dx);
  intcoords[2] = coordMin + IntType((realcoords[2] - xmin[2])/dx);
}

//------------------------------------------------------------------------------
// Convert int coordinates to reals.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
void
QuantizedTessellation3d<IntType, RealType>::
dequantize(const IntType* intcoords, RealType* realcoords) const {
  const RealType dx = length/(coordMax - coordMin);
  realcoords[0] = xmin[0] + (intcoords[0] - coordMin)*dx;
  realcoords[1] = xmin[1] + (intcoords[1] - coordMin)*dx;
  realcoords[2] = xmin[2] + (intcoords[2] - coordMin)*dx;
}

//------------------------------------------------------------------------------
// Read out the current quantized tessellation state to a regular tessellation.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
void
QuantizedTessellation3d<IntType, RealType>::
fillTessellation(Tessellation<3, RealType>& mesh) const {
  POLY_ASSERT2(false, "Implement me!");
}

//------------------------------------------------------------------------------
// Internal method to construct once we have the bounds set.
// Assumes xmin, xmax, and length are already set.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
void
QuantizedTessellation3d<IntType, RealType>::
construct(const std::vector<RealType>& points) {
  POLY_ASSERT(points.size() % 3 == 0);
  POLY_ASSERT(std::abs(xmax[0] - xmin[0] - length) < 1e-10*length);
  POLY_ASSERT(std::abs(xmax[1] - xmin[1] - length) < 1e-10*length);
  POLY_ASSERT(std::abs(xmax[2] - xmin[2] - length) < 1e-10*length);
  const int numGenerators = points.size()/3;
  generators.resize(numGenerators);
  for (unsigned i = 0; i < numGenerators; ++i) {
    this->quantize(&points[3*i], &generators[i].x);
    generators[i].index = i;
  }
  POLY_ASSERT2(false, "Implement 3D guard generators!");
}

//------------------------------------------------------------------------------
// Instantiate versions we know we need.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType> 
IntType
QuantizedTessellation3d<IntType, RealType>::coordMin = std::numeric_limits<IntType>::min()/3;

template<typename IntType, typename RealType> 
IntType
QuantizedTessellation3d<IntType, RealType>::coordMax = std::numeric_limits<IntType>::max()/3;

template struct QuantizedTessellation3d<int, double>;
template struct QuantizedTessellation3d<int64_t, double>;
}
