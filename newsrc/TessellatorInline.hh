// #include "makeBoxPLC.hh"
#include "findBoundaryElements.hh"
#include "SiloWriter.hh"
// #include "snapToBoundary.hh"
// #include "polytope_geometric_utilities.hh"

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
  auto& Q = Quantizer<nDim>::instance();
  if (!Q.m_init) {
    Q.init(points);
    m_init = true;
  }
  // Pre-conditions
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(points.size() > 0);
  POLY_ASSERT(points.size() % nDim == 0);

  // Invoke the descendant method to fill the quant mesh.
  QuantTessellation<nDim> quantmesh(points);
  this->tessellateQuantized(quantmesh);

  // Copy the QuantTessellation to the output.
  quantmesh.fillTessellation(mesh);

  // Fill in the boundary elements.
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
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
           const PLC<nDim>& geometry,
           Tessellation<nDim, RealType>& mesh) const {
  auto& Q = Quantizer<nDim>::instance();
  if (!Q.m_init) {
    Q.init(PLCpoints);
    m_init = true;
  }
  // Pre-conditions
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(points.size() > 0);
  POLY_ASSERT(points.size() % nDim == 0);

  // Invoke the descendant method to fill the quant mesh.
  QuantTessellation<nDim> quantmesh(points);
  QuantPLC<nDim> qplc(geometry, PLCpoints);
  // Remove any external points
  quantmesh.cullExternalPoints(qplc);
  this->tessellateQuantized(quantmesh);

  // Clip against the boundary.
  // Remove non-facet points and merge collinear facets
  quantmesh.clipTessellation(qplc, *this);

  // Copy the QuantTessellation to the output.
  quantmesh.fillTessellation(mesh);

  // Fill in the boundary elements.
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
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

}
