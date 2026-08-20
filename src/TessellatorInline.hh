#include "findBoundaryElements.hh"

namespace polytope {

//------------------------------------------------------------------------
// Tessellate (unbounded)
//------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
void
Tessellator<Dimension, RealType>::
tessellate(const std::vector<RealType>& points,
           Tessellation<Dimension, RealType>& mesh) {
  if (points.size() == 0) {
    return;
  }
  auto& Q = Quantizer<Dimension>::instance();
  if (!Q.m_init) {
    Q.init(points);
  }
  // Pre-conditions
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(points.size() > 0);
  POLY_ASSERT(points.size() % Dimension == 0);

  // Invoke the descendant method to fill the quant mesh.
  QuantTessellation<Dimension> quantmesh(points);
  this->tessellateQuantized(quantmesh);

  // Convert back to physical space.
  quantmesh.fillTessellation(mesh);

  // Fill in the boundary elements.
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
}

//------------------------------------------------------------------------
// Tessellate in a PLC.
//------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
void
Tessellator<Dimension, RealType>::
tessellate(const std::vector<RealType>& points,
           const std::vector<RealType>& PLCpoints,
           const PLC<Dimension>& geometry,
           Tessellation<Dimension, RealType>& mesh) {
  if (points.size() == 0) {
    return;
  }
  auto& Q = Quantizer<Dimension>::instance();
  if (!Q.m_init) {
    Q.init(PLCpoints);
  }
  // Pre-conditions
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(points.size() > 0);
  POLY_ASSERT(points.size() % Dimension == 0);

  // Invoke the descendant method to fill the quant mesh.
  QuantTessellation<Dimension> quantmesh(points);
  QuantPLC<Dimension> qplc(geometry, PLCpoints);
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

//------------------------------------------------------------------------
// Manual tessellation if only a single generator is provided.
//------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline void
Tessellator<Dimension, RealType>::
singleNodeTessellate(QuantTessellation<Dimension>& result) {
  if constexpr (Dimension == 2) {
    const auto& Q = Quantizer<2>::instance();
    result.cells.resize(1);

    // Map canonical edges to face indices for orientation tracking
    edge::EdgeToFaceMap edgeToFace;

    // Map QuantizedPoint coordinates to node indices for deduplication
    std::map<QuantizedPoint<2>, int> node2id;

    // Add nodes for the box extent and keep track of their indices
    auto cornerIndices = addBoxPoints(Q, node2id, result.nodes);
    const int N = 4;
    BoxSides side;
    for (int i = 0; i < N; ++i) {
      auto point0 = cornerIndices[side.corner(i)];
      auto point1 = cornerIndices[side.corner((i+1)%N)];
      int signedFaceIndex = edge::addOrientedEdge(point0, point1, result.faces, edgeToFace);
      result.cells[0].push_back(signedFaceIndex);
    }
  }
}

}
