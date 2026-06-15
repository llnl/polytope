//-----------------------------------------------------------------------------//
// QuantTessellation<3> specializations
//-----------------------------------------------------------------------------//

#include "QuantTessellation.hh"
#include "Intersections.hh"
#include "Tessellator.hh"
#include <map>
#include <set>

namespace polytope {

// template<>
// std::vector<std::vector<Point3<HashKey<3>::IntType>>>
// QuantTessellation<3>::getCell(const unsigned cellIndx) const {
//   auto numFaces = m_cells[cellIndx];
//   IntCell polygon;
//   // TODO: implement me
//   return polygon;
// }

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
  // Storage for clipped cells (will replace m_cells)
  std::vector<std::vector<int>> newCells;
  std::vector<IntPoint> newNodes;
  std::vector<std::vector<int>> newFaces;

  // Storage for generator points corresponding to surviving cells
  std::vector<IntPoint> newPoints;
  std::vector<CoordHash> newHashes;

  // Map from IntPoint to index in newNodes (for vertex deduplication)
  std::map<IntPoint, int> nodeMap;

  auto getOrCreateNode = [&](const IntPoint& p) -> int {
    auto it = nodeMap.find(p);
    if (it != nodeMap.end()) {
      return it->second;
    }
    int newIdx = newNodes.size();
    newNodes.push_back(p);
    nodeMap[p] = newIdx;
    return newIdx;
  };

  // Clip each cell against all PLC boundary planes
  for (size_t iCell = 0; iCell < m_cells.size(); ++iCell) {
    const auto& cellFaceIndices = m_cells[iCell];

    // Extract cell geometry: collect all vertices and faces for this cell
    std::set<int> cellNodeSet;
    std::vector<std::vector<int>> cellFaces;

    for (int faceIdx : cellFaceIndices) {
      // Handle signed face indices (negative = inverted orientation)
      int absFaceIdx = (faceIdx < 0) ? ~faceIdx : faceIdx;
      POLY_ASSERT(absFaceIdx < m_faces.size());

      const auto& faceNodes = m_faces[absFaceIdx];
      cellFaces.push_back(faceNodes);

      // Collect all nodes used by this cell
      for (int nodeIdx : faceNodes) {
        cellNodeSet.insert(nodeIdx);
      }
    }

    // Build local vertex array and remap face indices to local
    std::vector<int> globalToLocal(m_nodes.size(), -1);
    std::vector<IntPoint> localVertices;
    localVertices.reserve(cellNodeSet.size());

    for (int globalIdx : cellNodeSet) {
      int localIdx = localVertices.size();
      globalToLocal[globalIdx] = localIdx;
      localVertices.push_back(m_nodes[globalIdx]);
    }

    // Remap face indices to local vertex indices
    std::vector<std::vector<int>> localFaces;
    localFaces.reserve(cellFaces.size());
    for (const auto& face : cellFaces) {
      std::vector<int> localFace;
      localFace.reserve(face.size());
      for (int globalIdx : face) {
        POLY_ASSERT(globalToLocal[globalIdx] >= 0);
        localFace.push_back(globalToLocal[globalIdx]);
      }
      localFaces.push_back(localFace);
    }

    // Clip against each PLC boundary plane
    bool fullyClipped = false;
    for (size_t iFacet = 0; iFacet < QPLC.facets.size() && !fullyClipped; ++iFacet) {
      // Get plane definition from PLC facet
      const auto& planeNormal = QPLC.m_normals[iFacet];
      const auto& planePoint = QPLC.m_points[QPLC.facets[iFacet][0]];

      bool fullyRetained;
      clipPolyhedronByPlane(localVertices, localFaces,
                           planePoint, planeNormal,
                           fullyClipped, fullyRetained);
    }

    // If cell survived clipping, add it to the new tessellation
    if (!fullyClipped && !localVertices.empty() && !localFaces.empty()) {
      // Map clipped vertices to global node indices
      std::vector<int> localToNewGlobal(localVertices.size());
      for (size_t i = 0; i < localVertices.size(); ++i) {
        localToNewGlobal[i] = getOrCreateNode(localVertices[i]);
      }

      // Create new faces with global indices
      std::vector<int> newCellFaceIndices;
      for (const auto& localFace : localFaces) {
        std::vector<int> globalFace;
        globalFace.reserve(localFace.size());
        for (int localIdx : localFace) {
          POLY_ASSERT(localIdx < localToNewGlobal.size());
          globalFace.push_back(localToNewGlobal[localIdx]);
        }

        // Add face to global face list
        int newFaceIdx = newFaces.size();
        newFaces.push_back(globalFace);
        newCellFaceIndices.push_back(newFaceIdx);
      }

      // Add clipped cell
      newCells.push_back(newCellFaceIndices);

      // Keep the generator point for this surviving cell
      // This maintains the invariant: newCells[i] corresponds to newPoints[i]
      newPoints.push_back(m_points[iCell]);
      newHashes.push_back(m_hashes[iCell]);
    }
    // If fullyClipped: generator is discarded (not added to newPoints)
  }

  // Replace tessellation data with clipped version
  m_nodes = std::move(newNodes);
  m_faces = std::move(newFaces);
  m_cells = std::move(newCells);

  // Update generator points to match surviving cells
  m_points = std::move(newPoints);
  m_hashes = std::move(newHashes);

  // Update bounding box to match remaining generators
  if (!m_points.empty()) {
    m_loBounds = m_points[0];
    m_hiBounds = m_points[0];
    for (const auto& p : m_points) {
      m_loBounds = m_loBounds.minElements(p);
      m_hiBounds = m_hiBounds.maxElements(p);
    }
  } else {
    // All generators were clipped away
    m_loBounds = m_Q.maxCoord;
    m_hiBounds = -m_loBounds;
  }
}

//------------------------------------------------------------------------------
// Fill 3D tessellation mesh
//------------------------------------------------------------------------------
template<>
void
QuantTessellation<3>::fillTessellation(TessellationType& mesh) const {
  const unsigned numNodes = m_nodes.size();
  const unsigned numFaces = m_faces.size();
  const unsigned numCells = m_points.size();  // Number of generators

  // Allocate space for mesh data
  // In 3D: nodes are stored as [x0, y0, z0, x1, y1, z1, ...]
  mesh.nodes.resize(3 * numNodes);
  mesh.faces.resize(numFaces);
  mesh.faceCells.resize(numFaces);
  mesh.cells = m_cells;
  POLY_ASSERT2(m_cells.size() == numCells, "Differing number of cells and generator points");

  // Dequantize nodes from integer coordinates to real coordinates
  for (unsigned i = 0; i != numNodes; ++i) {
    RealPoint rp = m_Q.dequantize(m_nodes[i]);
    mesh.nodes[3*i]     = rp.x;
    mesh.nodes[3*i + 1] = rp.y;
    mesh.nodes[3*i + 2] = rp.z;
  }

  // Copy face topology (each face has 3+ nodes in 3D - triangular or polygonal)
  for (unsigned i = 0; i != numFaces; ++i) {
    mesh.faces[i].resize(m_faces[i].size());
    for (unsigned j = 0; j != m_faces[i].size(); ++j) {
      mesh.faces[i][j] = m_faces[i][j];
    }
  }

  // Build faceCells connectivity: for each cell, mark which faces touch it
  // Cells store signed face indices where negative means inverted orientation.
  for (unsigned i = 0; i != numCells; ++i) {
    const unsigned nf = mesh.cells[i].size();
    for (unsigned j = 0; j != nf; ++j) {
      auto k = mesh.cells[i][j];
      if (k < 0) {
        // Negative index: inverted face orientation
        POLY_ASSERT2(~k < numFaces, k << " " << ~k << " " << numFaces);
        mesh.faceCells[~k].push_back(~i);
      } else {
        // Positive index: normal face orientation
        POLY_ASSERT2(k < numFaces, k << " " << numFaces);
        mesh.faceCells[k].push_back(i);
      }
    }
  }
}

//------------------------------------------------------------------------------
// Explicit instantiation
//------------------------------------------------------------------------------
template class QuantTessellation<3>;

}
