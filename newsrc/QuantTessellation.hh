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
  using IntPoint = Point<Dimension, IntType>;
  using RealPoint = Point<Dimension, RealType>;
  using QuantizerType = Quantizer<Dimension>;
  using TessellationType = Tessellation<Dimension, RealType>;
  using IntCell = typename Cell<Dimension, IntType>::CellType;
  using RealCell = typename Cell<Dimension, RealType>::CellType;

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

  IntCell getCell (const unsigned cellIndx) const {
    return Cell<Dimension, IntType>::extractCell(m_nodes, m_cells[cellIndx], m_faces);
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
  // Dimension-specific methods (implemented in separate .cc files)
  //------------------------------------------------------------------------------

  // Dequantize and fill the tessellation type after tessellation
  void fillTessellation(TessellationType& mesh) const;

  // Clip tessellation against PLC boundary planes
  void clipTessellation(const QuantPLC<Dimension>& QPLC,
                        const Tessellator<Dimension, double>& tessellator);

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
