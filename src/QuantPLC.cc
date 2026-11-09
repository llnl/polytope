//-----------------------------------------------------------------------------//
// QuantPLC definitions
//-----------------------------------------------------------------------------//

#include <map>
#include <set>

#include "polytope_internal.hh"
#include "QuantPLC.hh"
#include "GeomUtils.hh"
#include "Intersections.hh"
#include "EdgeUtils.hh"
#include "Shapes.hh"
#ifdef POLYTOPE_ENABLE_QHULL
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#endif

namespace polytope {

template<int Dimension>
void
QuantPLC<Dimension>::init(const std::vector<QuantizedPoint<Dimension>>& points,
                          const PLC<Dimension>& plc) {
  facets = plc.facets;
  holes = plc.holes;
  init(points);
}

// Initialize with a vector of quantized points. Be sure to retain the indices
template<int Dimension>
void
QuantPLC<Dimension>::init(const std::vector<QuantizedPoint<Dimension>>& points) {
  const auto& Q = m_Q;
  m_loBounds = Q.maxBound;
  m_hiBounds = -m_loBounds;
  m_reduced = false;
  m_convex = false;
  for (const auto& ip : points) {
    m_loBounds = m_loBounds.minElements(ip);
    m_hiBounds = m_hiBounds.maxElements(ip);
  }
  if (isValid(points)) {
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
QuantPLC<Dimension>::reduce(std::vector<QuantizedPoint<Dimension>>& points) {
  if (!isValid(points)) return;
  // Collect only unique point indices used in the PLC facets
  std::set<unsigned> indices;
  for (const auto& f : facets) {
    std::copy(f.begin(), f.end(), std::inserter(indices, indices.end()));
  }
  for (const auto& hole : holes) {
    for (const auto& f : hole) {
      std::copy(f.begin(), f.end(), std::inserter(indices, indices.end()));
    }
  }

  const auto& Q = m_Q;
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
  std::vector<QuantizedPoint<Dimension>> reducedPoints;
  reducedPoints.reserve(N);
  for (auto h : newHashes) {
    auto point = Q.decode(h);
    point.index = reducedPoints.size();
    reducedPoints.push_back(point);
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

  points = std::move(reducedPoints);
  m_loBounds = Q.maxBound;
  m_hiBounds = -m_loBounds;
  for (const auto& point : points) {
    m_loBounds = m_loBounds.minElements(point);
    m_hiBounds = m_hiBounds.maxElements(point);
  }
  m_reduced = true;
}

//------------------------------------------------------------------------------
// makeConvex functions
//------------------------------------------------------------------------------
template<int Dimension>
void
QuantPLC<Dimension>::makeConvex(const std::vector<QuantizedPoint<Dimension>>& a_points) {
  m_convex = false;
  facets.clear();
  holes.clear();
  init(a_points);
  if constexpr (Dimension == 2) {
    makeConvex2D<Dimension>(a_points);
  } else if constexpr (Dimension == 3) {
    makeConvex3D<Dimension>(a_points);
  }
}

// 2D
//---------------------
template<int Dimension>
template<int D>
std::enable_if_t<D == 2, void>
QuantPLC<Dimension>::makeConvex2D(const std::vector<QuantizedPoint<2>>& a_points) {
  const unsigned n = a_points.size();
  if (n < 3) return;
  m_convex = true;
#ifdef POLYTOPE_ENABLE_QHULL
  std::vector<double> q_points;
  q_points.reserve(2*n);
  for (auto& p : a_points) {
    q_points.push_back(static_cast<double>(p.x));
    q_points.push_back(static_cast<double>(p.y));
  }
  orgQhull::Qhull qh("2D points", 2, static_cast<int>(n), q_points.data(), "");
  orgQhull::QhullFacet qfacet = qh.facetList().first();
  orgQhull::QhullVertex from;
  qfacet = qfacet.nextFacet2d(&from);
  const int edgeCount = qh.facetList().count();
  facets.clear();
  facets.resize(edgeCount);
  for (int i = 0; i < edgeCount; ++i) {
    orgQhull::QhullVertex to;
    qfacet = qfacet.nextFacet2d(&to);
    facets[i] = {static_cast<unsigned>(from.point().id()),
                 static_cast<unsigned>(to.point().id())};
    from = to;
  }
#else
  // Non-Qhull method is not as optimized
  // Unhash all points to integer coordinates and pair with original indices
  std::vector<std::pair<QuantizedPoint<Dimension>, unsigned>> sortedPoints;
  sortedPoints.reserve(n);
  for (unsigned i = 0; i < n; ++i) {
    sortedPoints.push_back(std::make_pair(a_points[i], i));
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
  orderFacets();
#endif
}

// 3D
//---------------------
template<int Dimension>
template<int D>
std::enable_if_t<D == 3, void>
QuantPLC<Dimension>::makeConvex3D(const std::vector<QuantizedPoint<3>>& a_points) {
#ifdef POLYTOPE_ENABLE_QHULL
  const unsigned n = a_points.size();
  if (n == 0) return;
  m_convex = true;

  std::vector<double> q_points;
  q_points.reserve(3*n);
  for (auto& p : a_points) {
    q_points.push_back(static_cast<double>(p.x));
    q_points.push_back(static_cast<double>(p.y));
    q_points.push_back(static_cast<double>(p.z));
  }
  orgQhull::Qhull qh("", 3, static_cast<int>(n), q_points.data(), "Qt");

  facets.clear();
  facets.resize(qh.facetList().size());
  auto i = 0;
  // Extract and reformat the Qhull points
  orgQhull::QhullFacetListIterator qfacets(qh.facetList());
  while (qfacets.hasNext()) {
    const orgQhull::QhullFacet qfacet = qfacets.next();
    orgQhull::QhullVertexSetIterator vertices(qfacet.vertices());
    std::array<orgQhull::QhullPoint, 3> p = {vertices.next().point(),
                                             vertices.next().point(),
                                             vertices.next().point()};
    // Ensure order is CCW using cross product and the plane normal
    const double ux = p[1][0] - p[0][0];
    const double uy = p[1][1] - p[0][1];
    const double uz = p[1][2] - p[0][2];
    const double vx = p[2][0] - p[0][0];
    const double vy = p[2][1] - p[0][1];
    const double vz = p[2][2] - p[0][2];
    const orgQhull::QhullHyperplane normal = qfacet.hyperplane();
    const double orientation =
      (uy * vz - uz * vy) * normal[0] +
      (uz * vx - ux * vz) * normal[1] +
      (ux * vy - uy * vx) * normal[2];
    std::array<unsigned, 3> pid = {static_cast<unsigned>(p[0].id()),
                                   static_cast<unsigned>(p[1].id()),
                                   static_cast<unsigned>(p[2].id())};
    if (orientation < 0.0) {
      std::swap(pid[1], pid[2]);
    }
    facets[i++] = {pid[0], pid[1], pid[2]};
  }
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
  orderEdgeLoop(facets);

  for (auto& hole : holes) {
    orderEdgeLoop(hole);
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
QuantPLC<Dimension>::within(const RealPoint& point,
                            const std::vector<QuantizedPoint<Dimension>>& points) const {
  const auto& Q = m_Q;
  if (!Q.inBounds(point)) {
    return false;
  }
  QuantizedPoint<Dimension> ip = Q.quantize(point);
  return within(ip, points);
}

template<int Dimension>
bool
QuantPLC<Dimension>::within(const QuantizedPoint<Dimension>& point,
                            const std::vector<QuantizedPoint<Dimension>>& points) const {
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
    return within2D<Dimension>(point, points);
  } else if constexpr (Dimension == 3) {
    return within3D<Dimension>(point, points);
  }
}

template<int Dimension>
template<int D>
std::enable_if_t<D == 2, bool>
QuantPLC<Dimension>::within2D(const QuantizedPoint<Dimension>& point,
                              const std::vector<QuantizedPoint<Dimension>>& points) const {
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
QuantPLC<Dimension>::within3D(const QuantizedPoint<Dimension>& /*point*/,
                              const std::vector<QuantizedPoint<Dimension>>& /*points*/) const {
  // Implement this
  return false;
}

//------------------------------------------------------------------------------
// Instantiate versions we know we need.
//------------------------------------------------------------------------------
template class QuantPLC<2>;
template class QuantPLC<3>;

}
