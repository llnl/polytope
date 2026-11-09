//-----------------------------------------------------------------------------//
// QuantPLC
//
// A Piecewise Linear Complex whose vertex coordinates are supplied by its
// caller. Facet indices refer to positions in that point vector. Call init()
// again if the coordinates change; this refreshes the cached bounds.
// Combines functionality from:
//   - PLC: Stores only points used in facets with index remapping
//   - convexHull_2d/3d: Can compute convex hulls using quantized coordinates
//   - Quantizer: Uses the selected key encoding for deduplication
//
// Functionality provided:
//   1. Automatic deduplication via hash-based comparison (no fuzzy tolerance)
//   2. Efficient spatial queries via encoded-key ordering
//   3. Unified interface for both hull construction and PLC reduction
//-----------------------------------------------------------------------------//

#ifndef __Polytope_QuantPLC__
#define __Polytope_QuantPLC__

#include <type_traits>
#include <utility>

#include "PLC.hh"
#include "QuantizedKeyTraits.hh"
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
  using Wide = QuantizedKey<Dimension>;
  using WidePoint = Point<Dimension, Wide>;
  using Quant = Quantizer<Dimension>;
  using QuantizedCell = Cell<Dimension, QuantizedCoordinate<Dimension>>;

  explicit QuantPLC(const Quant& Q): m_Q(Q) {}
  virtual ~QuantPLC() = default;

  QuantPLC(const QuantPLC&) = default;
  QuantPLC(QuantPLC&&) = default;

  // References cannot be rebound. Assignment therefore retains this object's
  // quantizer and is only valid between PLCs built with the same context.
  QuantPLC& operator=(const QuantPLC& other) {
    POLY_VERIFY2(&m_Q == &other.m_Q,
                 "Cannot assign QuantPLCs that use different quantizers");
    if (this != &other) {
      PLC<Dimension>::operator=(other);
      m_reduced = other.m_reduced;
      m_convex = other.m_convex;
      m_loBounds = other.m_loBounds;
      m_hiBounds = other.m_hiBounds;
    }
    return *this;
  }

  QuantPLC& operator=(QuantPLC&& other) {
    POLY_VERIFY2(&m_Q == &other.m_Q,
                 "Cannot assign QuantPLCs that use different quantizers");
    if (this != &other) {
      PLC<Dimension>::operator=(std::move(other));
      m_reduced = other.m_reduced;
      m_convex = other.m_convex;
      m_loBounds = other.m_loBounds;
      m_hiBounds = other.m_hiBounds;
    }
    return *this;
  }

  QuantPLC(const Quant& Q,
           const std::vector<QuantizedPoint<Dimension>>& points,
           const PLC<Dimension>& plc = PLC<Dimension>()):
    PLC<Dimension>(plc), m_Q(Q) { init(points); }

  // Quantize real coordinates into a caller-owned vector.
  static std::vector<QuantizedPoint<Dimension>>
  quantizePoints(const Quant& Q, const std::vector<RealPoint>& realPoints) {
    std::vector<QuantizedPoint<Dimension>> result;
    result.reserve(realPoints.size());
    for (const auto& point : realPoints) {
      auto quantized = Q.quantize(point);
      quantized.index = result.size();
      result.push_back(quantized);
    }
    return result;
  }

  void init(const std::vector<QuantizedPoint<Dimension>>& points,
            const PLC<Dimension>& plc);
  void init(const std::vector<QuantizedPoint<Dimension>>& points);

  // Reduce the caller-owned vector and remap all facet and hole indices.
  // Do not pass a generator vector with parallel hash/cell arrays here.
  void reduce(std::vector<QuantizedPoint<Dimension>>& points);

  void makeConvex(const std::vector<QuantizedPoint<Dimension>>& a_points);

  // Order facets to form proper boundaries
  void orderFacets();

  bool within(const QuantizedPoint<Dimension>& point,
              const std::vector<QuantizedPoint<Dimension>>& points) const;

  bool within(const RealPoint& point,
              const std::vector<QuantizedPoint<Dimension>>& points) const;

  QuantizedCell getCell(const std::vector<QuantizedPoint<Dimension>>& points) const {
    return QuantizedCell(points, facets);
  }

  std::vector<QuantizedCell>
  getHolePoints(const std::vector<QuantizedPoint<Dimension>>& points) const {
    std::vector<QuantizedCell> holePoints;
    holePoints.reserve(holes.size());
    for (const auto& hole : holes) {
      holePoints.push_back(QuantizedCell(points, hole));
    }
    return holePoints;
  }

  bool isValid(const std::vector<QuantizedPoint<Dimension>>& points) const {
    if (points.size() < Dimension + 1) return false;
    return true;
  }

  static bool convexPLCIntersection(const QuantPLC<Dimension>& a,
                                    const std::vector<QuantizedPoint<Dimension>>& aPoints,
                                    const QuantPLC<Dimension>& b,
                                    const std::vector<QuantizedPoint<Dimension>>& bPoints) {
    if (!a.isValid(aPoints) || !b.isValid(bPoints)) {
      return false;
    }
    POLY_ASSERT2((a.m_convex && b.m_convex), "Must call makeConvex() on both inputs");
    if constexpr (Dimension == 2) {
      return convexIntersect<QuantizedCoordinate<Dimension>>(a.getCell(aPoints), b.getCell(bPoints));
    } else if constexpr (Dimension == 3) {
      // Implement this
      return false;
    }
  }

  //------------------------------------------------------------------------
  //! Functions used for testing
  //------------------------------------------------------------------------
  bool equals(const std::vector<QuantizedPoint<Dimension>>& points,
              const QuantPLC<Dimension>& other,
              const std::vector<QuantizedPoint<Dimension>>& otherPoints) const {
    return getCell(points) == other.getCell(otherPoints);
  }

  //------------------------------------------------------------------------
  // Member data
  //------------------------------------------------------------------------
  using PLC<Dimension>::facets; // Facets as vertex index lists
  using PLC<Dimension>::holes; // Holes (each hole is a list of facets)

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

  const Quant& quantizer() const { return m_Q; }

private:
  const Quant& m_Q;
  template<int D = Dimension>
  std::enable_if_t<D == 2, void> makeConvex2D(const std::vector<QuantizedPoint<2>>& a_points);

  template<int D = Dimension>
  std::enable_if_t<D == 2, void> orderFacets2D();

  template<int D = Dimension>
  std::enable_if_t<D == 3, void> orderFacets3D();

  template<int D = Dimension>
  std::enable_if_t<D == 3, void> makeConvex3D(const std::vector<QuantizedPoint<3>>& a_points);

  template<int D = Dimension>
  std::enable_if_t<D == 2, bool>
  within2D(const QuantizedPoint<Dimension>& point,
           const std::vector<QuantizedPoint<Dimension>>& points) const;

  template<int D = Dimension>
  std::enable_if_t<D == 3, bool>
  within3D(const QuantizedPoint<Dimension>& point,
           const std::vector<QuantizedPoint<Dimension>>& points) const;
};

// Serialization
// This is only used for the distributed tessellator so only some things
// need to be serialized
template<int Dimension>
struct Serializer<QuantPLC<Dimension>> {
  static void serializeImpl(const QuantPLC<Dimension>& value,
                            std::vector<char>& buffer) {
    serialize(value.facets, buffer);
    serialize(value.m_convex, buffer);
    serialize(value.m_loBounds, buffer);
    serialize(value.m_hiBounds, buffer);
  }

  static void deserializeImpl(QuantPLC<Dimension>& value,
                              std::vector<char>::const_iterator& bufItr,
                              const std::vector<char>::const_iterator endItr) {
    deserialize(value.facets, bufItr, endItr);
    deserialize(value.m_convex, bufItr, endItr);
    deserialize(value.m_loBounds, bufItr, endItr);
    deserialize(value.m_hiBounds, bufItr, endItr);
  }
};
}
#endif
