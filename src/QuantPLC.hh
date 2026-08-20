//-----------------------------------------------------------------------------//
// QuantPLC
//
// A Piecewise Linear Complex that stores points in both quantized (integer)
// and real coordinate spaces, using Morton curve hashing for efficient
// spatial queries. Combines functionality from:
//   - PLC: Stores only points used in facets with index remapping
//   - convexHull_2d/3d: Can compute convex hulls using quantized coordinates
//   - MortonKeyTraits/Quantizer: Uses Morton keys for deduplication
//
// Functionality provided:
//   1. Automatic deduplication via hash-based comparison (no fuzzy tolerance)
//   2. Efficient spatial queries via Morton curve ordering
//   3. Unified interface for both hull construction and PLC reduction
//-----------------------------------------------------------------------------//

#ifndef __Polytope_QuantPLC__
#define __Polytope_QuantPLC__

#include <type_traits>

#include "PLC.hh"
#include "MortonKeyTraits.hh"
#include "Quantizer.hh"
#include "Intersections.hh"
#include "Cell.hh"
#include "Serializer.hh"

namespace polytope {

// TODO: Have this somehow inherit the routines in
// QuantTessellation since there is lots of overlap
template<int Dimension>
class QuantPLC : public PLC<Dimension> {
public:
  using RealType = double;
  using RealPoint = Point<Dimension, RealType>;
  using Wide = MortonKey<Dimension>;
  using WidePoint = Point<Dimension, Wide>;
  using Quant = Quantizer<Dimension>;
  using QuantizedCell = typename Cell<Dimension, QuantizedCoordinate<Dimension>>::CellType;

  QuantPLC() = default;
  virtual ~QuantPLC() = default;

  QuantPLC(const PLC<Dimension>& plc,
           const std::vector<RealType>& allpoints);

  explicit QuantPLC(const std::vector<RealType>& allpoints);

  QuantPLC(const PLC<Dimension>& plc,
           const std::vector<QuantizedPoint<Dimension>>& quantizedPoints);

  void init(const PLC<Dimension>& plc,
            const std::vector<RealType>& allpoints);

  void init(const std::vector<RealType>& allpoints);

  void init(const PLC<Dimension>& plc,
            const std::vector<QuantizedPoint<Dimension>>& quantizedPoints);

  void init(const std::vector<QuantizedPoint<Dimension>>& quantizedPoints);

  // Remove any degenerate points
  void removeDegeneracies();

  // Reduce to only the points used in the boundary facets, if they exist
  void reduce();

  void makeConvex();

  // Order facets to form proper boundaries
  void orderFacets();

  bool within(const QuantizedPoint<Dimension>& point) const;

  bool within(const RealPoint& point) const;

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
    const auto& Q = Quant::instance();
    std::vector<RealPoint> realPoints;
    realPoints.reserve(points.size());
    for (const auto& p : points) {
      realPoints.push_back(Q.dequantize(p));
    }
    return flattenCoords(realPoints);
  }

  QuantizedCell getCell() const {
    return Cell<Dimension, QuantizedCoordinate<Dimension>>::extractCell(points, facets);
  }

  std::vector<QuantizedCell> getHolePoints() const {
    std::vector<QuantizedCell> holePoints;
    holePoints.reserve(holes.size());
    for (const auto& hole : holes) {
      holePoints.push_back(Cell<Dimension, QuantizedCoordinate<Dimension>>::extractCell(points, hole));
    }
    return holePoints;
  }

  bool isValid() const {
    if (points.size() < Dimension + 1) return false;
    return true;
  }

  static bool convexPLCIntersection(const QuantPLC<Dimension>& a,
                                    const QuantPLC<Dimension>& b) {
    if (!a.isValid() || !b.isValid()) {
      return false;
    }
    POLY_ASSERT2((a.m_convex && b.m_convex), "Must call makeConvex() on both inputs");
    if constexpr (Dimension == 2) {
      return convexIntersect<QuantizedCoordinate<Dimension>>(a.getCell(), b.getCell());
    } else if constexpr (Dimension == 3) {
      // Implement this
      return false;
    }
  }

  //------------------------------------------------------------------------
  //! Functions used for testing
  //------------------------------------------------------------------------
  std::vector<MortonKey<Dimension>> sortedHashes() const {
    std::vector<MortonKey<Dimension>> sorted(hashes);
    std::sort(sorted.begin(), sorted.end());
    return sorted;
  }

  static bool compareHashes(const QuantPLC<Dimension>& lhs,
                            const QuantPLC<Dimension>& rhs) {
    return lhs.sortedHashes() == rhs.sortedHashes();
  }

  std::vector<std::set<MortonKey<Dimension>>> facetHashSet() const {
    std::vector<std::set<MortonKey<Dimension>>> facetSet;
    for (const auto& f : facets) {
      facetSet.push_back(std::set<MortonKey<Dimension>>());
      for (const auto& idx : f) {
        facetSet.back().insert(hashes[idx]);
      }
    }
    return facetSet;
  }

  static bool compareFacets(const QuantPLC<Dimension>& lhs,
                            const QuantPLC<Dimension>& rhs) {
    auto set1 = lhs.facetHashSet();
    auto set2 = rhs.facetHashSet();
    if (set1.size() != set2.size()) {
      return false;
    }
    for (const auto& s1 : set1) {
      bool found = false;
      for (const auto& s2 : set2) {
        if (s1 == s2) {
          found = true;
          break;
        }
      }
      if (!found) return false;
    }
    return true;
  }

  //------------------------------------------------------------------------------
  // Member data
  //------------------------------------------------------------------------------
  using PLC<Dimension>::facets; // Facets as vertex index lists
  using PLC<Dimension>::holes; // Holes (each hole is a list of facets)
  std::vector<MortonKey<Dimension>> hashes;
  std::vector<QuantizedPoint<Dimension>> points;

  // Precomputed geometric properties (3D only)
  // std::vector<QuantizedPoint<Dimension>> m_normals; // Normalized normals for each facet
  // std::vector<WidePoint> m_faceCentroidSums; // Unnormalized centroid sums
  // std::vector<Wide> m_faceCentroidCounts;    // Vertex counts for each facet
  // WidePoint m_polyCentroidSum;               // Unnormalized polyhedron centroid
  // Wide m_polyCentroidWeight;                 // Weight for polyhedron centroid

  bool m_reduced = false;
  bool m_convex = false;
  // Local lower and upper bounding box coordinates
  QuantizedPoint<Dimension> m_loBounds;
  QuantizedPoint<Dimension> m_hiBounds;

private:
  template<int D = Dimension>
  std::enable_if_t<D == 2, void> makeConvex2D();

  template<int D = Dimension>
  std::enable_if_t<D == 2, void> orderFacets2D();

  template<int D = Dimension>
  std::enable_if_t<D == 3, void> orderFacets3D();

  template<int D = Dimension>
  std::enable_if_t<D == 3, void> makeConvex3D();

  template<int D = Dimension>
  std::enable_if_t<D == 2, bool>
  within2D(const QuantizedPoint<Dimension>& point) const;

  template<int D = Dimension>
  std::enable_if_t<D == 3, bool>
  within3D(const QuantizedPoint<Dimension>& point) const;
};

// Serialization
// This is only used for the distributed tessellator so only some things
// need to be serialized
template<int Dimension>
struct Serializer<QuantPLC<Dimension>> {
  static void serializeImpl(const QuantPLC<Dimension>& value,
                            std::vector<char>& buffer) {
    serialize(value.facets, buffer);
    serialize(value.points, buffer);
    serialize(value.m_convex, buffer);
    serialize(value.m_loBounds, buffer);
    serialize(value.m_hiBounds, buffer);
  }

  static void deserializeImpl(QuantPLC<Dimension>& value,
                              std::vector<char>::const_iterator& bufItr,
                              const std::vector<char>::const_iterator endItr) {
    deserialize(value.facets, bufItr, endItr);
    deserialize(value.points, bufItr, endItr);
    deserialize(value.m_convex, bufItr, endItr);
    deserialize(value.m_loBounds, bufItr, endItr);
    deserialize(value.m_hiBounds, bufItr, endItr);
  }
};
}
#endif
