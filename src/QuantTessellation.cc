//-----------------------------------------------------------------------------//
// QuantTessellation definitions
//-----------------------------------------------------------------------------//

//#include "polytope_internal.hh"
#include "QuantTessellation.hh"
#include "Intersections.hh"
#include <map>
#include <set>

namespace polytope {

template<typename Dimension>
QuantTessellation<Dimension>::
QuantTessellation(const Quant& Q,
                  const std::vector<RealType>& genpoints) :
  m_Q(Q) {
  m_loBounds.one();
  m_hiBounds.zero();
  // Extract the unrolled coordinates
  std::vector<RealPoint> rpoints = extractCoords<Dimension, RealType>(genpoints);

  auto N = rpoints.size();
  m_hashes.reserve(N);
  m_points.reserve(N);
  for (const auto& rp : rpoints) {
    auto ip = m_Q.quantize(rp);
    m_loBounds = m_loBounds.minElements(ip);
    m_hiBounds = m_hiBounds.maxElements(ip);
    m_hashes.push_back(m_Q.hash(ip));
    m_points.push_back(ip);
  }
}

template<typename Dimension>
std::vector<RealPoint>
QuantTessellation<Dimension>::getRealPoints() const {
  std::vector<RealPoint> realPoints;
  realPoints.reserve(m_points.size());
  for (const auto& p : m_points) {
    realPoints.push_back(p.template type_cast<RealType>());
  }
  return realPoints;
}

template<typename Dimension>
std::vector<IntPoint>
QuantTessellation<Dimension>::getIntPoints() const {
  std::vector<IntPoint> intPoints(m_points);
  return intPoints;
}

template<typename Dimension>
void
QuantTessellation<Dimension>::
fillTessellation(TessellationType& mesh) const {
  if constexpr (Dimension == 2) {
    fillTessellation2D(mesh);
  } else if constexpr (Dimension == 3) {
    fillTessellation3D(mesh);
  }
}

template<int Dimension>
template<int D>
std::enable_if_t<D == 2, void>
QuantTessellation<Dimension>::
fillTessellation2D(TessellationType& mesh) const {
  const unsigned numNodes = m_nodes.size();
  const unsigned numFaces = m_faces.size();
  const unsigned numCells = m_points.size();  // Number of generators

  // Allocate space for mesh data
  // In 2D: nodes are stored as [x0, y0, x1, y1, ...]
  mesh.nodes.resize(2 * numNodes);
  mesh.faces.resize(numFaces, std::vector<unsigned>(2));
  mesh.faceCells.resize(numFaces);
  mesh.cells = m_cells;

  // Dequantize nodes from integer coordinates to real coordinates
  for (unsigned i = 0; i != numNodes; ++i) {
    RealPoint rp = m_Q.dequantize(m_nodes[i]);
    mesh.nodes[2*i]     = rp.x;
    mesh.nodes[2*i + 1] = rp.y;
  }

  // Copy face topology (each face has 2 nodes in 2D)
  for (unsigned i = 0; i != numFaces; ++i) {
    POLY_ASSERT(m_faces[i].size() == 2);
    mesh.faces[i][0] = m_faces[i][0];
    mesh.faces[i][1] = m_faces[i][1];
  }

  // Build faceCells connectivity: for each cell, mark which faces touch it
  // In the old implementation, cells stored signed face indices where
  // negative meant inverted orientation. The new implementation stores
  // positive indices, so we need to adapt if needed.
  for (unsigned i = 0; i != numCells; ++i) {
    const unsigned nf = mesh.cells[i].size();
    for (unsigned j = 0; j != nf; ++j) {
      int k = mesh.cells[i][j];
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

template<int Dimension>
template<int D>
std::enable_if_t<D == 3, void>
QuantTessellation<Dimension>::
fillTessellation3D(TessellationType& mesh) const {
  // TODO: Implement me
}

//------------------------------------------------------------------------------
// Explicit instantiations
//------------------------------------------------------------------------------
template class QuantTessellation<2>;
template class QuantTessellation<3>;

}
