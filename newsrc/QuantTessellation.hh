//-----------------------------------------------------------------------------//
// QuantTessellation
//
//-----------------------------------------------------------------------------//

#ifndef __Polytope_QuantTessellation__
#define __Polytope_QuantTessellation__

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

// Forward declaration
template<int Dimension, typename RealType>
class Tessellator;

// TODO: Make this inherit from Tessellation class
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
  QuantTessellation(const std::vector<RealType>& genpoints) {
    init(genpoints);
  }

  // Construct a smaller instance, useful during clipping
  QuantTessellation(const std::vector<IntPoint>& qgenpoints) {
    init(qgenpoints);
  }

  void init(const std::vector<RealType>& genpoints) {
    const auto& Q = QuantizerType::instance();
    m_loBounds = Q.maxCoord;
    m_hiBounds = -m_loBounds;
    // Extract the unrolled coordinates
    auto rpoints = extractCoords<Dimension, RealType>(genpoints);

    auto N = rpoints.size();
    hashes.reserve(N);
    points.reserve(N);
    size_t i = 0;
    for (const auto& rp : rpoints) {
      auto ip = Q.quantize(rp);
      POLY_ASSERT2(Q.inQBounds(ip), "Point provided that exceeds quantizer bounds");
      ip.index = i++;
      m_loBounds = m_loBounds.minElements(ip);
      m_hiBounds = m_hiBounds.maxElements(ip);
      hashes.push_back(Q.hash(ip));
      points.push_back(ip);
    }
    sortByHash();
  }

  void init(const std::vector<IntPoint>& qgenpoints) {
    points = qgenpoints;
    const auto& Q = QuantizerType::instance();
    m_loBounds = Q.maxCoord;
    m_hiBounds = -m_loBounds;
    auto N = points.size();
    hashes.reserve(N);
    unsigned i = 0;
    for (auto& ip : points) {
      ip.index = i++;
      m_loBounds = m_loBounds.minElements(ip);
      m_hiBounds = m_hiBounds.maxElements(ip);
      hashes.push_back(Q.hash(ip));
    }
  }

  void clear() {
    points.clear();
    hashes.clear();
    nodes.clear();
    faces.clear();
    cells.clear();
  }

  //------------------------------------------------------------------------------
  // Common methods
  //------------------------------------------------------------------------------

  // Returns quantized points cast as reals to give to the tessellator
  std::vector<RealPoint> getRealPoints() const {
    std::vector<RealPoint> realPoints;
    realPoints.reserve(points.size());
    for (const auto& p : points) {
      realPoints.push_back(p.template type_cast<RealType>());
    }
    return realPoints;
  }

  // Returns dequantized points cast as a flattened vector of reals
  std::vector<RealType> getRealCoords() const {
    const auto& Q = QuantizerType::instance();
    std::vector<RealPoint> realPoints;
    realPoints.reserve(points.size());
    for (const auto& p : points) {
      realPoints.push_back(Q.dequantize(p));
    }
    return flattenCoords(realPoints);
  }

  // Returns quantized points
  const std::vector<IntPoint> getIntPoints() const { return points; }

  IntCell getCell(const unsigned cellIndex) const {
    return Cell<Dimension, IntType>::extractCell(nodes, cells[cellIndex], faces);
  }

  // Reorder points and cells based on sorted hashes for deterministic output
  void sortByHash() {
    const auto numPoints = points.size();
    if (numPoints == 0) return;

    // Create index vector and sort by hash
    std::vector<unsigned> sortedIndices(numPoints);
    std::iota(sortedIndices.begin(), sortedIndices.end(), 0);
    std::sort(sortedIndices.begin(), sortedIndices.end(),
              [this](unsigned a, unsigned b) { return hashes[a] < hashes[b]; });

    // Create mapping from old index to new index
    std::vector<unsigned> oldToNew(numPoints);
    for (unsigned i = 0; i < numPoints; ++i) {
      oldToNew[sortedIndices[i]] = i;
    }

    // Reorder points and hashes
    std::vector<IntPoint> newPoints(numPoints);
    std::vector<CoordHash> newHashes(numPoints);
    for (unsigned i = 0; i < numPoints; ++i) {
      newPoints[i] = points[sortedIndices[i]];
      newPoints[i].index = i;
      newHashes[i] = hashes[sortedIndices[i]];
    }
    points = std::move(newPoints);
    hashes = std::move(newHashes);

    // Reorder cells array using the same permutation
    if (cells.size() > 0) {
      std::vector<std::vector<int>> newCells(numPoints);
      for (unsigned i = 0; i < numPoints; ++i) {
        newCells[i] = std::move(cells[sortedIndices[i]]);
      }
      cells = std::move(newCells);
    }
  }

  //------------------------------------------------------------------------------
  // Common methods for all dimensions
  //------------------------------------------------------------------------------

  // Remove unused nodes and faces (not referenced by any cells)
  void compactUnusedNodesAndFaces() {
    const unsigned numCells = points.size();

    // Step 1: Find all faces referenced by cells
    std::set<unsigned> usedFaces;
    for (unsigned i = 0; i < numCells; ++i) {
      for (auto faceIdx : cells[i]) {
        // Handle signed face indices (negative means inverted orientation)
        unsigned absFaceIdx = (faceIdx < 0) ? ~faceIdx : faceIdx;
        usedFaces.insert(absFaceIdx);
      }
    }

    // Step 2: Find all nodes referenced by used faces
    std::set<unsigned> usedNodes;
    for (auto faceIdx : usedFaces) {
      for (auto nodeIdx : faces[faceIdx]) {
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
      IntPoint node = nodes[oldIdx];
      node.index = oldToNewNode[oldIdx];
      newNodes.push_back(node);
    }

    std::vector<std::vector<unsigned>> newFaces;
    newFaces.reserve(usedFaces.size());
    for (auto oldIdx : usedFaces) {
      std::vector<unsigned> face;
      face.reserve(faces[oldIdx].size());
      for (auto nodeIdx : faces[oldIdx]) {
        face.push_back(oldToNewNode[nodeIdx]);
      }
      newFaces.push_back(face);
    }

    // Step 5: Remap cell face indices
    for (unsigned i = 0; i < numCells; ++i) {
      for (auto& faceIdx : cells[i]) {
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
    nodes = std::move(newNodes);
    faces = std::move(newFaces);
  }

  //------------------------------------------------------------------------------
  // Dimension-specific methods (implemented in separate .cc files)
  //------------------------------------------------------------------------------

  // Dequantize and fill the tessellation type after tessellation
  void fillTessellation(TessellationType& mesh);

  // Clip tessellation against PLC boundary planes
  void clipTessellation(const QuantPLC<Dimension>& QPLC,
                        Tessellator<Dimension, double>& tessellator);

  // Remove any external generator points
  void cullExternalPoints(const QuantPLC<Dimension>& QPLC);

  // Compare cells with convex hull
  bool cellIntersectsHull(const QuantPLC<Dimension>& QPLC,
                          const unsigned cellIndex) const;

  //------------------------------------------------------------------------------
  // Output method
  //------------------------------------------------------------------------------
  friend std::ostream& operator<<(std::ostream& s, const QuantTessellation& mesh) {
    for (int i = 0; i < mesh.cells.size(); ++i) {
      s << mesh.getCell(i);
    }
    return s;
  }

  //------------------------------------------------------------------------------
  // Member data
  //------------------------------------------------------------------------------
  // Generator points
  std::vector<CoordHash> hashes;
  std::vector<IntPoint> points;
  // Local lower and upper bounding box coordinates
  IntPoint m_loBounds;
  IntPoint m_hiBounds;
  std::vector<IntPoint> nodes; // Nodes that make up the Voronoi
  std::vector<std::vector<unsigned>> faces; // Faces made up of indices into nodes
  std::vector<std::vector<int>> cells; // Cells made up of indices into faces
};

} // namespace polytope
#endif
