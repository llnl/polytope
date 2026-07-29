//-----------------------------------------------------------------------------//
// QuantTessellation<3> specializations
//-----------------------------------------------------------------------------//

#include "QuantTessellation.hh"
#include "Intersections.hh"
#include "Tessellator.hh"
#include <map>
#include <set>

namespace polytope {

// Check if any cells intersect a convex hull
template<>
bool
QuantTessellation<3>::cellIntersectsHull(const QuantPLC<3>& QPLC,
                                         const unsigned cellIndex) const {
  auto plc_cell = QPLC.getFacetPoints();
  auto qcell = getCell(cellIndex);
  // TODO: Implement me
  //return convexIntersection(plc_cell, qcell);
  return true;
}

// Remove any generator points that are outside our clipping region
template<>
void
QuantTessellation<3>::cullExternalPoints(const QuantPLC<3>& QPLC) {
  POLY_CONTRACT_VAR(QPLC);
  // TODO: Implement me
}

//------------------------------------------------------------------------------
// Clip 3D tessellation against PLC boundary planes
//
// Clips each convex cell (polyhedron) in the tessellation against all
// boundary planes defined by the QuantPLC. Cells that are fully clipped away
// are removed, and their corresponding generator points are also removed.
// The resulting tessellation contains only the portion inside the PLC boundary
// with proper cell-generator correspondence maintained.
//
// Algorithm:
//   1. For each cell, extract its local geometry (vertices + faces)
//   2. Clip the cell against each PLC boundary plane sequentially
//   3. If the cell survives, remap its vertices to a new global vertex list
//   4. Keep the generator point for surviving cells, discard for clipped cells
//   5. Build new global face and cell lists from the clipped geometry
//------------------------------------------------------------------------------
template<>
void
QuantTessellation<3>::clipTessellation(const QuantPLC<3>& QPLC,
                                       const Tessellator<3, double>& tessellator) {
  POLY_CONTRACT_VAR(QPLC);
  POLY_CONTRACT_VAR(tessellator);
  // 3D has not been implemented
}

//------------------------------------------------------------------------------
// Fill 3D tessellation mesh
//------------------------------------------------------------------------------
template<>
void
QuantTessellation<3>::fillTessellation(TessellationType& mesh) {
  auto& Q = Quantizer<3>::instance();
  compactUnusedNodesAndFaces();
  const unsigned numNodes = m_nodes.size();
  const unsigned numFaces = m_faces.size();
  const unsigned numCells = m_points.size();  // Number of generators

  // Allocate space for mesh data
  // In 3D: nodes are stored as [x0, y0, z0, x1, y1, z1, ...]
  mesh.nodes.resize(numNodes);
  mesh.faces.resize(numFaces);
  mesh.cells = m_cells;
  POLY_ASSERT2(m_cells.size() == numCells, "Differing number of cells and generator points");

  for (unsigned i = 0; i < numCells; ++i) {
    RealPoint rp = Q.dequantize(m_points[i]);
    mesh.points[i].x = rp.x;
    mesh.points[i].y = rp.y;
    mesh.points[i].z = rp.z;
  }

  // Dequantize nodes from integer coordinates to real coordinates
  for (unsigned i = 0; i < numNodes; ++i) {
    RealPoint rp = Q.dequantize(m_nodes[i]);
    mesh.nodes[i].x = rp.x;
    mesh.nodes[i].y = rp.y;
    mesh.nodes[i].z = rp.z;
  }

  // Copy face topology (each face has 3+ nodes in 3D - triangular or polygonal)
  for (unsigned i = 0; i < numFaces; ++i) {
    mesh.faces[i].resize(m_faces[i].size());
    for (unsigned j = 0; j < m_faces[i].size(); ++j) {
      mesh.faces[i][j] = m_faces[i][j];
    }
  }
  mesh.computeFaceCells();
}

//------------------------------------------------------------------------------
// Explicit instantiation
//------------------------------------------------------------------------------
template class QuantTessellation<3>;

}
