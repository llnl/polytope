//-----------------------------------------------------------------------------//
// QuantTessellation
//
//-----------------------------------------------------------------------------//

#ifndef __Polytope_QuantTessellation__
#define __Polytope_QuantTessellation__

#include <vector>
#include <algorithm>
#include <numeric>
#include <string>
#include "Quantizer.hh"
#include "Tessellation.hh"
#include "QuantPLC.hh"
#include "Shapes.hh"
#include "Communicator.hh"
#include "Intersections.hh"

namespace polytope {

// Forward declaration
template<int Dimension, typename RealType>
class Tessellator;

// TODO: Make this inherit from Tessellation class
template<int Dimension>
class QuantTessellation : public Tessellation<Dimension, typename HashKey<Dimension>::IntType> {
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
  virtual ~QuantTessellation() {};

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

  QuantTessellation(const std::vector<std::vector<CoordHash>>& rankHashes) {
    init(rankHashes);
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

  // Initialize or extend the generator points
  void init(const std::vector<std::vector<CoordHash>>& rankHashes) {
    faces.clear();
    nodes.clear();
    cells.clear();
    faceCells.clear();
    const auto& Q = Quantizer<Dimension>::instance();
    auto nranks = rankHashes.size();
    auto Ntotal = points.size();
    for (auto& v : rankHashes) {
      Ntotal += v.size();
    }
    // This implies the cell ranks were not filled initially
    // Fill them now
    if (cellRank.size() != points.size()) {
      cellRank.assign(points.size(), Communicator::getRank());
    }
    points.reserve(Ntotal);
    hashes.reserve(Ntotal);
    cellRank.reserve(Ntotal);
    unsigned i = Ntotal;
    for (auto source = 0u; source < nranks; ++source) {
      for (auto& ch : rankHashes[source]) {
        auto ip = Q.unhash(ch);
        ip.index = i++;
        hashes.push_back(ch);
        points.push_back(ip);
        cellRank.push_back(source);
      }
    }
    sortByHash();
  }

  void clear() {
    Tessellation<Dimension, IntType>::clear();
    hashes.clear();
  }

  //------------------------------------------------------------------------------
  // Common methods
  //------------------------------------------------------------------------------

  // Returns quantized points cast as doubles to give to the tessellator
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

  virtual IntCell getCell(const unsigned cellIndex) const override {
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
    if (cellRank.size() > 0) {
      std::vector<int> newCellRank(numPoints);
      for (unsigned i = 0; i < numPoints; ++i) {
        newCellRank[i] = cellRank[sortedIndices[i]];
      }
      cellRank = std::move(newCellRank);
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

  void makeConvexHull() {
    if (convexHull.m_convex) {
      return;
    }
    PLC<Dimension> emptyPLC;
    // Check if hull will be valid
    convexHull.init(emptyPLC, points);
    convexHull.makeConvex();
  }

  //------------------------------------------------------------------------------
  // Parallel methods
  //------------------------------------------------------------------------------

  // Extract the visible generators
  std::vector<CoordHash> visibleGenerators();

  // Determine which ranks generators are neighbors
  std::set<int> neighboringRanks() {
    this->computeFaceCells();
    int rank = Communicator::getRank();
    std::set<int> neighbors;
    for (const auto& fc : faceCells) {
      if (fc.size() < 2) continue;
      auto af0 = fc[0] < 0 ? ~fc[0] : fc[0];
      auto af1 = fc[1] < 0 ? ~fc[1] : fc[1];
      const auto& ranki = cellRank[af0];
      const auto& rankj = cellRank[af1];
      if (ranki == rank && rankj != rank) neighbors.insert(rankj);
      if (rankj == rank && ranki != rank) neighbors.insert(ranki);
    }
    return neighbors;
  }

  void filterToLocalGenerators() {
    const auto N = points.size();
    auto rank = Communicator::getRank();
    std::vector<IntPoint> newPoints;
    std::vector<CoordHash> newHashes;
    std::vector<std::vector<int>> newCells;
    newPoints.reserve(N);
    newHashes.reserve(N);
    newCells.reserve(N);
    for (auto i = 0; i < N; ++i) {
      if (cellRank[i] == rank) {
        newPoints.push_back(points[i]);
        newHashes.push_back(hashes[i]);
        newCells.push_back(cells[i]);
      }
    }
    points = std::move(newPoints);
    hashes = std::move(newHashes);
    cells = std::move(newCells);
    compactUnusedNodesAndFaces();
  }

  //------------------------------------------------------------------------------
  // Used for debugging purposes. Creates a file of relevant generator points
  //------------------------------------------------------------------------------

  // Eject only certain generators
  void ejectEscapePod(std::string filename,
                      const std::vector<unsigned>& genPoints,
                      const QuantPLC<Dimension>& QPLC,
                      const std::string& tessellatorName = "");
  // Eject all generators
  void ejectEscapePod(std::string filename,
                      const QuantPLC<Dimension>& QPLC,
                      const std::string& tessellatorName = "") {
    std::vector<unsigned> genPoints(points.size());
    for (auto i = 0u; i < points.size(); ++i) {
      genPoints[i] = i;
    }
    ejectEscapePod(filename, genPoints, QPLC, tessellatorName);
  }

  void loadEscapePod(std::string filename,
                     QuantPLC<Dimension>& QPLC);
  void loadEscapePod(std::string filename,
                     QuantPLC<Dimension>& QPLC,
                     std::string& tessellatorName);

  bool cellIntersectsHull(const QuantPLC<Dimension>& QPLC,
                          const unsigned index) {
    auto qcell = getCell(index);
    auto plc_cell = QPLC.getCell();
    return convexBoundaryIntersect<IntType>(qcell, plc_cell);
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
  bool m_isEscapePod = false;
  std::vector<CoordHash> hashes;
  // Local lower and upper bounding box coordinates
  IntPoint m_loBounds;
  IntPoint m_hiBounds;
  using Tessellation<Dimension, IntType>::points; // Generator points
  using Tessellation<Dimension, IntType>::nodes; // Nodes that make up the Voronoi
  using Tessellation<Dimension, IntType>::faces; // Faces made up of indices into nodes
  using Tessellation<Dimension, IntType>::cells; // Cells made up of indices into faces
  using Tessellation<Dimension, IntType>::faceCells;
  using Tessellation<Dimension, IntType>::cellRank;
  QuantPLC<Dimension> convexHull;
};

} // namespace polytope
#endif
