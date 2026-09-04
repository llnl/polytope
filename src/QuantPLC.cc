//-----------------------------------------------------------------------------//
// QuantPLC definitions
//-----------------------------------------------------------------------------//

#include "polytope_internal.hh"
#include "QuantPLC.hh"
#include "GeomUtils.hh"
#include "Intersections.hh"
#ifdef POLYTOPE_ENABLE_QHULL
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#endif
#include <map>
#include <set>
#include "Shapes.hh"

namespace polytope {

template<int Dimension>
QuantPLC<Dimension>::QuantPLC(const PLC<Dimension>& plc,
                              const std::vector<RealType>& allpoints) :
  PLC<Dimension>(plc) {
  init(allpoints);
}

template<int Dimension>
QuantPLC<Dimension>::QuantPLC(const PLC<Dimension>& plc,
                              const std::vector<QuantizedPoint<Dimension>>& quantizedPoints) :
  PLC<Dimension>(plc) {
  init(quantizedPoints);
}

template<int Dimension>
QuantPLC<Dimension>::
QuantPLC(const std::vector<RealType>& allpoints) :
  QuantPLC(PLC<Dimension>(), allpoints) { }

template<int Dimension>
void
QuantPLC<Dimension>::init(const PLC<Dimension>& plc,
                          const std::vector<RealType>& allpoints) {
  facets = plc.facets;
  holes = plc.holes;
  init(allpoints);
}

template<int Dimension>
void
QuantPLC<Dimension>::init(const std::vector<RealType>& allpoints) {
  const auto& Q = Quantizer<Dimension>::instance();
  m_loBounds = Q.maxBound;
  m_hiBounds = -m_loBounds;

  // Extract the unrolled coordinates
  std::vector<RealPoint> rpoints = extractCoords<Dimension, RealType>(allpoints);

  auto N = rpoints.size();
  points.reserve(N);
  size_t i = 0;
  for (const auto& rp : rpoints) {
    auto ip = Q.quantize(rp);
    ip.index = i++;
    m_loBounds = m_loBounds.minElements(ip);
    m_hiBounds = m_hiBounds.maxElements(ip);
    points.push_back(ip);
  }
  if (isValid()) {
    POLY_ASSERT2(m_loBounds < m_hiBounds,
                 "Provided coplanar or collinear or degenerate points to the QuantPLC");
    orderFacets();
  }
}

template<int Dimension>
void
QuantPLC<Dimension>::init(const PLC<Dimension>& plc,
                          const std::vector<QuantizedPoint<Dimension>>& quantizedPoints) {
  facets = plc.facets;
  holes = plc.holes;
  init(quantizedPoints);
}

template<int Dimension>
void
QuantPLC<Dimension>::init(const std::vector<QuantizedPoint<Dimension>>& quantizedPoints) {
  const auto& Q = Quantizer<Dimension>::instance();
  m_loBounds = Q.maxBound;
  m_hiBounds = -m_loBounds;

  auto N = quantizedPoints.size();
  points.reserve(N);
  size_t i = 0;
  for (const auto& ip : quantizedPoints) {
    auto nip(ip);
    nip.index = i++;
    m_loBounds = m_loBounds.minElements(nip);
    m_hiBounds = m_hiBounds.maxElements(nip);
    points.push_back(nip);
  }
  if (isValid()) {
    POLY_ASSERT2(m_loBounds < m_hiBounds,
                 "Provided coplanar or collinear or degenerate points to the QuantPLC");
    orderFacets();
  }
}

//------------------------------------------------------------------------------
// Remove any points not associated with a facet and order the facets.
//------------------------------------------------------------------------------
template<int Dimension>
void
QuantPLC<Dimension>::reduce() {
  // Collect only unique point indices used in the PLC facets
  std::set<unsigned> indices;
  for (const auto& f : facets) {
    std::copy(f.begin(), f.end(), std::inserter(indices, indices.end()));
  }

  const auto& Q = Quantizer<Dimension>::instance();
  unsigned newIdx = 0;
  std::map<unsigned, unsigned> old2new;
  std::map<QuantizedKey<Dimension>, unsigned> hashToIndex;
  std::vector<QuantizedKey<Dimension>> newHashes;
  for(int oldIdx : indices) {
    auto hash = Q.encode(points[oldIdx]);
    if (hashToIndex.find(hash) != hashToIndex.end()) {
      old2new[oldIdx] = hashToIndex[hash];
    } else {
      old2new[oldIdx] = newIdx;
      hashToIndex[hash] = newIdx++;
      newHashes.push_back(hash);
    }
  }

  auto N = indices.size();
  points.clear();
  points.reserve(N);
  for (auto h : newHashes) {
    points.push_back(Q.decode(h));
  }

  // Remap facet indices
  for (auto& f : facets) {
    for (auto& idx : f) {
      idx = old2new[idx];
    }
  }
  for (auto& hole : holes) {
    for (auto& f : hole) {
      for (auto& idx : f) {
        idx = old2new[idx];
      }
    }
  }

  orderFacets();
  m_reduced = true;
}

//------------------------------------------------------------------------------
// makeConvex functions
//------------------------------------------------------------------------------
template<int Dimension>
void
QuantPLC<Dimension>::makeConvex() {
  if constexpr (Dimension == 2) {
    makeConvex2D<Dimension>();
  } else if constexpr (Dimension == 3) {
    makeConvex3D<Dimension>();
  }
}

// 2D
//---------------------
template<int Dimension>
template<int D>
std::enable_if_t<D == 2, void>
QuantPLC<Dimension>::makeConvex2D() {
  const unsigned n = points.size();
  if (n == 0) return;
  m_convex = true;
  // TODO: Replace this with QHull call instead

  // Unhash all points to integer coordinates and pair with original indices
  std::vector<std::pair<QuantizedPoint<Dimension>, unsigned>> sortedPoints;
  sortedPoints.reserve(n);
  for (unsigned i = 0; i < n; ++i) {
    sortedPoints.push_back(std::make_pair(points[i], i));
  }

  // Sort by x-coordinate, then y-coordinate (using Point's operator<)
  std::sort(sortedPoints.begin(), sortedPoints.end(),
            [](const std::pair<QuantizedPoint<Dimension>, unsigned>& a,
               const std::pair<QuantizedPoint<Dimension>, unsigned>& b) {
              return a.first < b.first;
            });

  // Check if points are collinear
  bool collinear = true;
  if (n > 2) {
    for (unsigned i = 2; i < n && collinear; ++i) {
      collinear = (aboveBelow(sortedPoints[0].first,
                              sortedPoints[1].first,
                              sortedPoints[i].first) == 0);
    }
  }

  // If collinear, the hull is just a line segment from first to last
  if (collinear) {
    facets.resize(1, std::vector<unsigned>(2));
    facets[0][0] = sortedPoints.front().second;
    facets[0][1] = sortedPoints.back().second;
    return;
  }

  // Non-collinear case: use monotone chain algorithm
  std::vector<unsigned> result(2 * n);
  unsigned k = 0;

  // Build lower hull
  for (unsigned i = 0; i < n; ++i) {
    while (k >= 2 &&
           aboveBelow(sortedPoints[result[k - 2]].first,
                      sortedPoints[result[k - 1]].first,
                      sortedPoints[i].first) >= 0) {
      k--;
    }
    result[k++] = i;
  }

  // Build upper hull
  unsigned t = k + 1;
  for (int i = n - 2; i >= 0; --i) {
    while (k >= t &&
           aboveBelow(sortedPoints[result[k - 2]].first,
                      sortedPoints[result[k - 1]].first,
                      sortedPoints[i].first) >= 0) {
      k--;
    }
    result[k++] = i;
  }

  POLY_ASSERT(k >= 4);
  POLY_ASSERT(result[0] == result[k - 1]);

  // Build facets from hull vertices (exclude duplicate last vertex)
  // The result array has k elements with result[0] == result[k-1]
  facets.clear();
  facets.reserve(k - 1);
  for (unsigned i = 0; i < k - 1; ++i) {
    facets.push_back(std::vector<unsigned>(2));
    facets.back()[0] = sortedPoints[result[i]].second;
    facets.back()[1] = sortedPoints[result[i + 1]].second;
  }
  if (isValid()) {
    // Remove any points that are not part of the convex hull
    reduce();
  }
}

// 3D
//---------------------
template<int Dimension>
template<int D>
std::enable_if_t<D == 3, void>
QuantPLC<Dimension>::makeConvex3D() {
#ifdef POLYTOPE_ENABLE_QHULL
  const unsigned n = points.size();
  if (n == 0) return;
  m_convex = true;

  std::vector<double> q_points = flattenCoords<3, double>(getRealQPoints());
  orgQhull::Qhull qh;
  qh.runQhull("", Dimension, n, q_points.data(), "Qt");

  facets.clear();
  facets.reserve(qh.facetList().size());
  // Extract and reformat the Qhull points
  for (auto f = qh.facetList().begin(); f != qh.facetList().end(); ++f) {
    auto vertices = f->vertices();
    facets.push_back(std::vector<unsigned>(vertices.size()));
    unsigned fIndx = 0;
    for (auto v = vertices.begin(); v != vertices.end(); ++v, fIndx++) {
      auto pid = (*v).point().id();
      facets.back()[fIndx] = pid;
    }
  }
  // Remove any points that are not part of the convex hull
  reduce();
#else
  if (Communicator::getRank() == Communicator::getRoot()) {
    std::cerr << "Must enable QHull to use 3D functionality" << std::endl;
  }
  Communicator::abort();
#endif
}

//------------------------------------------------------------------------------
// Order facets to form proper boundaries
//------------------------------------------------------------------------------
template<int Dimension>
void
QuantPLC<Dimension>::orderFacets() {
  if constexpr (Dimension == 2) {
    orderFacets2D();
  } else if constexpr (Dimension == 3) {
    orderFacets3D();
  }
}

//------------------------------------------------------------------------------
// Order 2D facets (edges) to form a closed loop
// Ensures facets[i][1] connects to facets[i+1][0]
// Also removes collinear points and joins adjacent collinear facets
//------------------------------------------------------------------------------
template<int Dimension>
template<int D>
std::enable_if_t<D == 2, void>
QuantPLC<Dimension>::orderFacets2D() {
  edge::orderEdgeLoop(facets);

  for (auto& hole : holes) {
    edge::orderEdgeLoop(hole);
    std::reverse(hole.begin(), hole.end());
  }
}

//------------------------------------------------------------------------------
// Order 3D facets to have consistent winding
// Ensures all facet normals point outward from the centroid
// Merges coplanar adjacent faces
// Works for arbitrary polygonal facets, not just triangles
//------------------------------------------------------------------------------
template<int Dimension>
template<int D>
std::enable_if_t<D == 3, void>
QuantPLC<Dimension>::orderFacets3D() {
  // Implement this
}

//------------------------------------------------------------------------------
// Intersection tests
//------------------------------------------------------------------------------
template<int Dimension>
bool
QuantPLC<Dimension>::within(const RealPoint& point) const {
  const auto& Q = Quantizer<Dimension>::instance();
  if (!Q.inBounds(point)) {
    return false;
  }
  QuantizedPoint<Dimension> ip = Q.quantize(point);
  return within(ip);
}

template<int Dimension>
bool
QuantPLC<Dimension>::within(const QuantizedPoint<Dimension>& point) const {
  POLY_ASSERT2(facets.size() > 0,
               "Must have facets before checking within");
  // First checks the bounds
  if (point < m_loBounds || point > m_hiBounds) {
    return false;
  }
  // Check if point is coincident with any vertices
  for (const auto& p : points) {
    if (p == point) {
      return true;
    }
  }
  if constexpr (Dimension == 2) {
    return within2D<Dimension>(point);
  } else if constexpr (Dimension == 3) {
    return within3D<Dimension>(point);
  }
}

template<int Dimension>
template<int D>
std::enable_if_t<D == 2, bool>
QuantPLC<Dimension>::within2D(const QuantizedPoint<Dimension>& point) const {
  // Build the outer boundary polygon from ordered facets (edges)
  // Since facets are now ordered, simply collect first vertex of each edge
  std::vector<QuantizedPoint<Dimension>> polygon;
  polygon.reserve(facets.size());
  for (const auto& f : facets) {
    polygon.push_back(points[f[0]]);
  }

  // Check if point is in or on the outer boundary
  if (polygon.size() < 3) return false;  // Degenerate case

  bool inOuter = pointInPolygon(point, polygon) || pointOnPolygon(point, polygon);
  if (!inOuter) return false;

  // Check if point is in any holes (if so, it's outside)
  for (const auto& hole : holes) {
    std::vector<QuantizedPoint<Dimension>> holePoly;
    holePoly.reserve(hole.size());
    for (const auto& f : hole) {
      holePoly.push_back(points[f[0]]);
    }
    // Points on the hole boundary do not count as inside the hole
    if (holePoly.size() >= 3 && pointInPolygon(point, holePoly)
        && !pointOnPolygon(point, holePoly)) {
      return false;
    }
  }

  return true;
}

template<int Dimension>
template<int D>
std::enable_if_t<D == 3, bool>
QuantPLC<Dimension>::within3D(const QuantizedPoint<Dimension>& /*point*/) const {
  // Implement this
  return false;
}

//------------------------------------------------------------------------------
// Instantiate versions we know we need.
//------------------------------------------------------------------------------
template class QuantPLC<2>;
template class QuantPLC<3>;

}
