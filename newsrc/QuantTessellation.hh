//-----------------------------------------------------------------------------//
// QuantTessellation
//
//-----------------------------------------------------------------------------//

#ifndef POLYTOPE_QUANTTESSELLATION_HH
#define POLYTOPE_QUANTTESSELLATION_HH

#include <vector>
#include <algorithm>
#include <numeric>
#include "HashKey.hh"
#include "Point.hh"
#include "Quantizer.hh"
#include "Tessellation.hh"
#include "QuantPLC.hh"
#include "Cell.hh"
#include "Shapes.hh"

namespace polytope {

template<int Dimension, typename RealType>
class Tessellator;

template<int Dimension>
class QuantTessellation {
public:
  using RealType = double;
  using CoordHash = typename HashKey<Dimension>::CoordHash;
  using IntType = typename HashKey<Dimension>::IntType;
  using HashType = typename HashKey<Dimension>::HashType;
  using IntPoint = Point<Dimension, IntType>;
  using RealPoint = Point<Dimension, RealType>;
  using QuantizerType = Quantizer<Dimension>;
  using TessellationType = Tessellation<Dimension, RealType>;
  using IntCell = typename Cell<Dimension, IntType>::CellType;
  using RealCell = typename Cell<Dimension, RealType>::CellType;

  QuantTessellation() = default;
  QuantTessellation& operator=(const QuantTessellation& other) = default;
  QuantTessellation(const QuantTessellation& other) = default;

  //------------------------------------------------------------------------------
  // Constructor - common for all dimensions
  //------------------------------------------------------------------------------
  QuantTessellation(const QuantizerType& Q,
                    const std::vector<RealType>& genpoints) :
    m_Q(Q) {
    m_loBounds = m_Q.maxCoord;
    m_hiBounds = -m_loBounds;
    // Extract the unrolled coordinates
    auto rpoints = extractCoords<Dimension, RealType>(genpoints);

    auto N = rpoints.size();
    m_hashes.reserve(N);
    m_points.reserve(N);
    size_t i = 0;
    for (const auto& rp : rpoints) {
      auto ip = m_Q.quantize(rp);
      POLY_ASSERT2(m_Q.inQBounds(ip), "Point provided that exceeds quantizer bounds");
      ip.index = i++;
      m_loBounds = m_loBounds.minElements(ip);
      m_hiBounds = m_hiBounds.maxElements(ip);
      m_hashes.push_back(m_Q.hash(ip));
      m_points.push_back(ip);
    }
    sortByHash();
  }

  // Construct a smaller instance, useful during clipping
  QuantTessellation(const std::vector<IntPoint>& qgenpoints,
                    const QuantTessellation& QT) :
    m_Q(QT.m_Q),
    m_points(qgenpoints),
    m_loBounds(QT.m_loBounds),
    m_hiBounds(QT.m_hiBounds) {
    auto N = qgenpoints.size();
    m_hashes.reserve(N);
    unsigned i = 0;
    for (auto& ip : m_points) {
      ip.index = i++;
      m_hashes.push_back(m_Q.hash(ip));
    }
  }

  void clear() {
    m_points.clear();
    m_hashes.clear();
    m_nodes.clear();
    m_faces.clear();
    m_cells.clear();
  }

  //------------------------------------------------------------------------------
  // Common methods
  //------------------------------------------------------------------------------

  // Returns quantized points cast as reals to give to the tessellator
  std::vector<RealPoint> getRealPoints() const {
    std::vector<RealPoint> realPoints;
    realPoints.reserve(m_points.size());
    for (const auto& p : m_points) {
      realPoints.push_back(p.template type_cast<RealType>());
    }
    return realPoints;
  }

  // Returns dequantized points cast as a flattened vector of reals
  std::vector<RealType> getRealCoords() const {
    std::vector<RealPoint> realPoints;
    realPoints.reserve(m_points.size());
    for (const auto& p : m_points) {
      realPoints.push_back(m_Q.dequantize(p));
    }
    return flattenCoords(realPoints);
  }

  // Returns quantized points
  const std::vector<IntPoint> getIntPoints() const { return m_points; }

  // Create guard generators
  std::vector<IntPoint> generateGuards() const {
    std::vector<IntPoint> guards;
    const auto len = m_Q.maxCoord - 1;
    IntPoint min;
    min.zero();
    return shapes::createBoxPoints(min, len);
  }

  // Return the quantized points with the guard generators appended at the end
  void guardGenerators(std::vector<IntPoint>& points) {
    auto guards = generateGuards();
    unsigned kk = m_points.size();
    for (auto& guard : guards) {
      guard.index = kk++;
      points.push_back(guard);
    }
  }

  IntCell getCell(const unsigned cellIndex) const {
    return Cell<Dimension, IntType>::extractCell(m_nodes, m_cells[cellIndex], m_faces);
  }

  // Reorder points and cells based on sorted hashes for deterministic output
  void sortByHash() {
    const auto numPoints = m_points.size();
    if (numPoints == 0) return;

    // Create index vector and sort by hash
    std::vector<unsigned> sortedIndices(numPoints);
    std::iota(sortedIndices.begin(), sortedIndices.end(), 0);
    std::sort(sortedIndices.begin(), sortedIndices.end(),
              [this](unsigned a, unsigned b) { return m_hashes[a] < m_hashes[b]; });

    // Create mapping from old index to new index
    std::vector<unsigned> oldToNew(numPoints);
    for (unsigned i = 0; i < numPoints; ++i) {
      oldToNew[sortedIndices[i]] = i;
    }

    // Reorder points and hashes
    std::vector<IntPoint> newPoints(numPoints);
    std::vector<CoordHash> newHashes(numPoints);
    for (unsigned i = 0; i < numPoints; ++i) {
      newPoints[i] = m_points[sortedIndices[i]];
      newPoints[i].index = i;
      newHashes[i] = m_hashes[sortedIndices[i]];
    }
    m_points = std::move(newPoints);
    m_hashes = std::move(newHashes);

    // Reorder cells array using the same permutation
    if (m_cells.size() > 0) {
      std::vector<std::vector<int>> newCells(numPoints);
      for (unsigned i = 0; i < numPoints; ++i) {
        newCells[i] = std::move(m_cells[sortedIndices[i]]);
      }
      m_cells = std::move(newCells);
    }
  }

  //------------------------------------------------------------------------------
  // Common methods for all dimensions
  //------------------------------------------------------------------------------

  // Remove unused nodes and faces (not referenced by any cells)
  void compactUnusedNodesAndFaces() {
    const unsigned numCells = m_points.size();

    // Step 1: Find all faces referenced by cells
    std::set<unsigned> usedFaces;
    for (unsigned i = 0; i < numCells; ++i) {
      for (auto faceIdx : m_cells[i]) {
        // Handle signed face indices (negative means inverted orientation)
        unsigned absFaceIdx = (faceIdx < 0) ? ~faceIdx : faceIdx;
        usedFaces.insert(absFaceIdx);
      }
    }

    // Step 2: Find all nodes referenced by used faces
    std::set<unsigned> usedNodes;
    for (auto faceIdx : usedFaces) {
      for (auto nodeIdx : m_faces[faceIdx]) {
        usedNodes.insert(nodeIdx);
      }
    }

    // Step 3: Create mapping from old indices to new compact indices
    std::map<unsigned, unsigned> oldToNewNode;
    unsigned newNodeIdx = 0;
    for (auto oldIdx : usedNodes) {
      oldToNewNode[oldIdx] = newNodeIdx++;
    }

    std::map<unsigned, unsigned> oldToNewFace;
    unsigned newFaceIdx = 0;
    for (auto oldIdx : usedFaces) {
      oldToNewFace[oldIdx] = newFaceIdx++;
    }

    // Step 4: Create compacted node and face arrays
    std::vector<IntPoint> newNodes;
    newNodes.reserve(usedNodes.size());
    for (auto oldIdx : usedNodes) {
      IntPoint node = m_nodes[oldIdx];
      node.index = oldToNewNode[oldIdx];
      newNodes.push_back(node);
    }

    std::vector<std::vector<int>> newFaces;
    newFaces.reserve(usedFaces.size());
    for (auto oldIdx : usedFaces) {
      std::vector<int> face;
      face.reserve(m_faces[oldIdx].size());
      for (auto nodeIdx : m_faces[oldIdx]) {
        face.push_back(oldToNewNode[nodeIdx]);
      }
      newFaces.push_back(face);
    }

    // Step 5: Remap cell face indices
    for (unsigned i = 0; i < numCells; ++i) {
      for (auto& faceIdx : m_cells[i]) {
        if (faceIdx < 0) {
          // Negative index: inverted face orientation
          unsigned oldFaceIdx = ~faceIdx;
          faceIdx = ~oldToNewFace[oldFaceIdx];
        } else {
          // Positive index: normal face orientation
          faceIdx = oldToNewFace[faceIdx];
        }
      }
    }

    // Step 6: Replace with compacted arrays
    m_nodes = std::move(newNodes);
    m_faces = std::move(newFaces);
  }

  //------------------------------------------------------------------------------
  // Dimension-specific methods (implemented in separate .cc files)
  //------------------------------------------------------------------------------

  // Dequantize and fill the tessellation type after tessellation
  void fillTessellation(TessellationType& mesh);

  // Clip tessellation against PLC boundary planes
  void clipTessellation(const QuantPLC<Dimension>& QPLC,
                        const Tessellator<Dimension, double>& tessellator);

  // Remove any external generator points
  void cullExternalPoints(const QuantPLC<Dimension>& QPLC);

  //------------------------------------------------------------------------------
  // Member data
  //------------------------------------------------------------------------------
  QuantizerType m_Q;
  // Generator points
  std::vector<CoordHash> m_hashes;
  std::vector<IntPoint> m_points;
  // Local lower and upper bounding box coordinates
  IntPoint m_loBounds;
  IntPoint m_hiBounds;
  std::vector<IntPoint> m_nodes; // Nodes that make up the Voronoi
  std::vector<std::vector<int>> m_faces; // Faces made up of indices into m_nodes
  std::vector<std::vector<int>> m_cells; // Cells made up of indices into m_faces
};

} // namespace polytope
#endif
