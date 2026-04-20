#include "clipQuantizedTessellation.hh"
#include "makeBoxPLC.hh"
#include "findBoundaryElements.hh"
#include "snapToBoundary.hh"

namespace polytope {

//----------------------------------------------------------------------------
// Tessellate (unbounded)
//------------------------------------------------------------------------------
template<int nDim, typename RealType>
inline
void
Tessellator<nDim, RealType>::
tessellate(const std::vector<RealType>& points,
           Tessellation<nDim, RealType>& mesh) const {

  // Pre-conditions
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(points.size() > 0);
  POLY_ASSERT(points.size() % nDim == 0);

  // Invoke the descendant method to fill the quant mesh.
  QuantizedTessellation quantmesh(points, points);
  this->tessellateQuantized(quantmesh);

  // Copy the QuantTessellation to the output.
  quantmesh.fillTessellation(mesh);

  // Fill in the boundary elements.
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
}

//----------------------------------------------------------------------------
// Tessellate in a rectangle.
//------------------------------------------------------------------------------
template<int nDim, typename RealType>
inline
void
Tessellator<nDim, RealType>::
tessellate(const std::vector<RealType>& points,
           RealType* low,
           RealType* high,
           Tessellation<nDim, RealType>& mesh) const {
  ReducedPLC<nDim, RealType> geometry;
  internal::makeBoxPLC(geometry, low, high);
  this->tessellate(points, geometry, mesh);
}

//----------------------------------------------------------------------------
// Tessellate in a PLC.
//------------------------------------------------------------------------------
template<int nDim, typename RealType>
inline
void
Tessellator<nDim, RealType>::
tessellate(const std::vector<RealType>& points,
           const std::vector<RealType>& PLCpoints,
           const PLC<nDim, RealType>& geometry,
           Tessellation<nDim, RealType>& mesh) const {
  
  // Pre-conditions
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(points.size() > 0);
  POLY_ASSERT(points.size() % nDim == 0);

  // Invoke the descendant method to fill the quant mesh.
  QuantizedTessellation quantmesh(points, PLCpoints);
  this->tessellateQuantized(quantmesh);

  // Clip against the boundary.
  clipQuantizedTessellation(quantmesh, PLCpoints, geometry, *this);

  // Copy the QuantTessellation to the output.
  quantmesh.fillTessellation(mesh);

  // Fill in the boundary elements.
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);

  // Snap exactly to the bounding PLC.
  snapToBoundary(mesh, PLCpoints, geometry, this->degeneracy());
}

//----------------------------------------------------------------------------
// Tessellate in a ReducedPLC.
//------------------------------------------------------------------------------
template<int nDim, typename RealType>
inline
void
Tessellator<nDim, RealType>::
tessellate(const std::vector<RealType>& points,
           const ReducedPLC<nDim, RealType>& geometry,
           Tessellation<nDim, RealType>& mesh) const {
  this->tessellate(points, geometry.points, geometry, mesh);
}

//------------------------------------------------------------------------------
// Do a tessellation with potential degenerate generators.
// Unbounded case.
//------------------------------------------------------------------------------
template<int nDim, typename RealType>
inline
std::vector<unsigned>
Tessellator<nDim, RealType>::
tessellateDegenerate(const std::vector<RealType>& points,
                     const RealType tol,
                     Tessellation<nDim, RealType>& mesh) const {

  // Compute the bounding box.
  RealType xmin[nDim], xmax[nDim];
  geometry::computeBoundingBox<nDim, RealType>(points, true, xmin, xmax);

  // Hash the input generators to a unique set.
  std::vector<RealType> uniquePoints;
  std::vector<unsigned> result;
  geometry::uniquePoints<nDim, RealType>(points, xmin, xmax, tol, uniquePoints, result);

  // Tessellate the unique set, and return the mapping.
  this->tessellate(uniquePoints, mesh);
  return result;
}

//------------------------------------------------------------------------------
// Do a tessellation with potential degenerate generators.
// Bounded by a box.
//------------------------------------------------------------------------------
template<int nDim, typename RealType>
inline
std::vector<unsigned>
Tessellator<nDim, RealType>::
tessellateDegenerate(const std::vector<RealType>& points,
                     RealType* low,
                     RealType* high,
                     const RealType tol,
                     Tessellation<nDim, RealType>& mesh) const {

  // Hash the input generators to a unique set.
  std::vector<RealType> uniquePoints;
  std::vector<unsigned> result;
  geometry::uniquePoints<nDim, RealType>(points, low, high, tol, uniquePoints, result);

  // Tessellate the unique set, and return the mapping.
  this->tessellate(uniquePoints, low, high, mesh);
  return result;
}

//------------------------------------------------------------------------------
// Do a tessellation with potential degenerate generators.
// PLC case.
//------------------------------------------------------------------------------
template<int nDim, typename RealType>
inline
std::vector<unsigned>
Tessellator<nDim, RealType>::
tessellateDegenerate(const std::vector<RealType>& points,
                     const std::vector<RealType>& PLCpoints,
                     const PLC<nDim, RealType>& geometry,
                     const RealType tol,
                     Tessellation<nDim, RealType>& mesh) const {

  // Compute the bounding box.
  RealType xmin[nDim], xmax[nDim];
  geometry::computeBoundingBox<nDim, RealType>(PLCpoints, true, xmin, xmax);

  // Hash the input generators to a unique set.
  std::vector<RealType> uniquePoints;
  std::vector<unsigned> result;
  geometry::uniquePoints<nDim, RealType>(points, xmin, xmax, tol, uniquePoints, result);

  // Tessellate the unique set, and return the mapping.
  this->tessellate(uniquePoints, PLCpoints, geometry, mesh);
  return result;
}

//------------------------------------------------------------------------------
// Do a tessellation with potential degenerate generators.
// ReducedPLC case.
//------------------------------------------------------------------------------
template<int nDim, typename RealType>
inline
std::vector<unsigned>
Tessellator<nDim, RealType>::
tessellateDegenerate(const std::vector<RealType>& points,
                     const ReducedPLC<nDim, RealType>& geometry,
                     const RealType tol,
                     Tessellation<nDim, RealType>& mesh) const {
  return this->tessellateDegenerate(points, geometry.points, geometry, tol, mesh);
}

//------------------------------------------------------------------------------
// Normalize points into a unit cube, optionally computing the bounds first.
//------------------------------------------------------------------------------
template<int nDim, typename RealType>
inline
std::vector<RealType>
Tessellator<nDim, RealType>::
computeNormalizedPoints(const std::vector<RealType>& points,
                        const std::vector<RealType>& PLCpoints,
                        const bool computeBounds,
                        RealType* low,
                        RealType* high) const {
  POLY_ASSERT(points.size() > 0);
  POLY_ASSERT(points.size() % nDim == 0);
  POLY_ASSERT(PLCpoints.empty() || PLCpoints.size() % nDim == 0);

  if (computeBounds) {
    geometry::computeBoundingBox<nDim, RealType>(points, false, low, high);
    if (!PLCpoints.empty()) {
      RealType plcLow[nDim], plcHigh[nDim];
      geometry::computeBoundingBox<nDim, RealType>(PLCpoints, false, plcLow, plcHigh);
      for (unsigned j = 0; j != nDim; ++j) {
        low[j] = std::min(low[j], plcLow[j]);
        high[j] = std::max(high[j], plcHigh[j]);
      }
    }
  }

  RealType length = RealType(0);
  for (unsigned j = 0; j != nDim; ++j) length = std::max(length, high[j] - low[j]);

  std::vector<RealType> result(points.size(), RealType(0));
  if (length > RealType(0)) {
    for (unsigned i = 0; i != points.size()/nDim; ++i) {
      for (unsigned j = 0; j != nDim; ++j) {
        result[nDim*i + j] = (points[nDim*i + j] - low[j])/length;
      }
    }
  }

  return result;
}

}
