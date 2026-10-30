#include "findBoundaryElements.hh"
#include "EdgeUtils.hh"

namespace polytope {

//------------------------------------------------------------------------
// Tessellate (unbounded)
//------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
Tessellation<Dimension, RealType>
Tessellator<Dimension, RealType>::
tessellate(const std::vector<Point<Dimension, RealType>>& points) {
  if (points.size() == 0) {
    return Tessellation<Dimension, RealType>();
  }
  initQ(points);

  // Pre-conditions
  POLY_ASSERT(points.size() > 0);

  // Invoke the descendant method to fill the quant mesh.
  QuantTessellation<Dimension> quantmesh(m_Q, points);
  quantmesh.tessellatorName = this->name();
  this->tessellateQuantized(quantmesh);

  // Convert back to physical space.
  Tessellation<Dimension, RealType> mesh;
  quantmesh.fillTessellation(mesh);

  // Fill in the boundary elements.
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  return mesh;
}

//------------------------------------------------------------------------
// Tessellate in a PLC.
//------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
Tessellation<Dimension, RealType>
Tessellator<Dimension, RealType>::
tessellate(const std::vector<Point<Dimension, RealType>>& points,
           const std::vector<Point<Dimension, RealType>>& PLCpoints,
           const PLC<Dimension>& geometry) {
  if (points.size() == 0) {
    return Tessellation<Dimension, RealType>();
  }
  initQ(points, PLCpoints);

  // Pre-conditions
  POLY_ASSERT(points.size() > 0);

  // Invoke the descendant method to fill the quant mesh.
  QuantTessellation<Dimension> quantmesh(m_Q, points);
  quantmesh.tessellatorName = this->name();
  QuantPLC<Dimension> qplc(m_Q, PLCpoints, geometry);
  // Remove any external points
  quantmesh.cullExternalPoints(qplc);
  this->tessellateQuantized(quantmesh);

  // Clip against the boundary.
  // Remove non-facet points and merge collinear facets
  quantmesh.clipTessellation(qplc, *this);

  // Copy the QuantTessellation to the output.
  Tessellation<Dimension, RealType> mesh;
  quantmesh.fillTessellation(mesh);

  // Fill in the boundary elements.
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  return mesh;
}

//------------------------------------------------------------------------
// Manual tessellation if only a single generator is provided.
//------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline void
Tessellator<Dimension, RealType>::
singleNodeTessellate(QuantTessellation<Dimension>& result) {
  if constexpr (Dimension == 2) {
    result.cells.resize(1);

    // Map canonical edges to face indices for orientation tracking
    EdgeToFaceMap edgeToFace;

    // Map QuantizedPoint coordinates to node indices for deduplication
    std::map<QuantizedPoint<2>, int> node2id;

    // Add nodes for the box extent and keep track of their indices
    auto cornerIndices = addBoxPoints(node2id, result.nodes, result.quantizer());
    const int N = 4;
    BoxSides side;
    for (int i = 0; i < N; ++i) {
      auto point0 = cornerIndices[side.corner(i)];
      auto point1 = cornerIndices[side.corner((i+1)%N)];
      int signedFaceIndex = addOrientedEdge(point0, point1, result.faces, edgeToFace);
      result.cells[0].push_back(signedFaceIndex);
    }
  }
}

}
