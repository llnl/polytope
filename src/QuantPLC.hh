//-----------------------------------------------------------------------------//
// QuantPLC
//
// A Piecewise Linear Complex that stores points in both quantized (integer)
// and real coordinate spaces, using Morton curve hashing for efficient
// spatial queries. Combines functionality from:
//   - ReducedPLC: Stores only points used in facets with index remapping
//   - convexHull_2d/3d: Can compute convex hulls using quantized coordinates
//   - HashKey/Quantizer: Uses Morton curve hashing for deduplication
//
// Key improvements over separate ReducedPLC + convexHull approach:
//   1. Automatic deduplication via hash-based comparison (no fuzzy tolerance)
//   2. Efficient spatial queries via Morton curve ordering
//   3. Unified interface for both hull construction and PLC reduction
//-----------------------------------------------------------------------------//

#ifndef POLYTOPE_QUANTPLC_HH
#define POLYTOPE_QUANTPLC_HH

#include <type_traits>

#include "PLC.hh"
#include "HashKey.hh"
#include "Quantizer.hh"

namespace polytope {

template<int Dimension>
class QuantPLC : public PLC<Dimension> {
public:
  using RealType = double;
  using CoordHash = typename HashKey<Dimension>::CoordHash;
  using IntType = typename HashKey<Dimension>::IntType;
  using IntPoint = Point<Dimension, IntType>;
  using RealPoint = Point<Dimension, RealType>;
  using Quant = Quantizer<Dimension>;
  using Hasher = HashKey<Dimension>;

  virtual ~QuantPLC() = default;

  QuantPLC(const PLC<Dimension>& plc,
           const Quant& Q,
           const std::vector<RealType>& allpoints);

  QuantPLC(const Quant& Q,
           const std::vector<RealType>& allpoints);

  // Remove any degenerate points
  void removeDegeneracies();

  // Reduce to only the points used in the boundary facets, if they exist
  void reduce();

  void makeConvex();

  // Order facets to form proper boundaries
  void orderFacets();

  bool within(const IntPoint& point) const;

  bool within(const RealPoint& point) const;

  //------------------------------------------------------------------------
  //! Functions used for testing
  //------------------------------------------------------------------------
  std::vector<CoordHash> sortedHashes() const {
    std::vector<CoordHash> sorted(m_hashes);
    std::sort(sorted.begin(), sorted.end());
    return sorted;
  }

  static bool compareHashes(const QuantPLC<Dimension>& lhs,
                            const QuantPLC<Dimension>& rhs) {
    return lhs.sortedHashes() == rhs.sortedHashes();
  }

  std::vector<std::set<CoordHash>> facetHashSet() const {
    std::vector<std::set<CoordHash>> facetSet;
    for (const auto& f : facets) {
      facetSet.push_back(std::set<CoordHash>());
      for (const auto& idx : f) {
        facetSet.back().insert(m_hashes[idx]);
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
  using PLC<Dimension>::facets;  // Facets as vertex index lists
  using PLC<Dimension>::holes;  // Holes (each hole is a list of facets)
  const Quant m_Q;
  std::vector<CoordHash> m_hashes;
  std::vector<IntPoint> m_points;
  bool m_reduced = false;
  bool m_convex_hull = false;
  // Local lower and upper bounding box coordinates
  IntPoint m_loBounds;
  IntPoint m_hiBounds;

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
  within2D(const IntPoint& point) const;

  template<int D = Dimension>
  std::enable_if_t<D == 3, bool>
  within3D(const IntPoint& point) const;
};

}
#endif
