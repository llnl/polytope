//-----------------------------------------------------------------------------//
// QuantPLC definitions
//-----------------------------------------------------------------------------//

#include "polytope_internal.hh"
#include "QuantPLC.hh"
#include "Intersections.hh"
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include <map>
#include <set>

namespace polytope {
namespace { // Anonymous namespace for internal helpers

//------------------------------------------------------------------------------
// sign of the Z coordinate of cross product : (p2 - p1)x(p3 - p1).
// Works directly with quantized integer coordinates.
//------------------------------------------------------------------------------
template<int Dimension, typename IntType>
int zcross_sign(const Point<Dimension, IntType>& p1,
                const Point<Dimension, IntType>& p2,
                const Point<Dimension, IntType>& p3) {
  // Use double precision to avoid overflow in the cross product calculation
  const double ztest =
    (double(p2.x) - double(p1.x))*(double(p3.y) - double(p1.y)) -
    (double(p2.y) - double(p1.y))*(double(p3.x) - double(p1.x));
  return (ztest < 0.0 ? -1 :
          ztest > 0.0 ?  1 :
                         0);
}

} // end anonymous namespace

template<int Dimension>
QuantPLC<Dimension>::QuantPLC(const PLC<Dimension>& plc,
                              const Quant& Q,
                              const std::vector<RealType>& allpoints) :
  PLC<Dimension>(plc),
  m_Q(Q) {
  m_loBounds.one();
  m_hiBounds.zero();

  // Extract the unrolled coordinates
  std::vector<RealPoint> rpoints = extractCoords<Dimension, RealType>(allpoints);

  auto N = rpoints.size();
  m_hashes.reserve(N);
  for (const auto& rp : rpoints) {
    auto ip = m_Q.quantize(rp);
    m_loBounds = m_loBounds.minElements(ip);
    m_hiBounds = m_hiBounds.maxElements(ip);
  }
  POLY_ASSERT2(m_loBounds < m_hiBounds,
               "Provided coplanar or collinear or degenerate points to the QuantPLC");
  removeDegeneracies();
}

template<int Dimension>
QuantPLC<Dimension>::
QuantPLC(const Quant& Q,
         const std::vector<RealType>& allpoints) :
  QuantPLC(PLC<Dimension>(), Q, allpoints) { }

//------------------------------------------------------------------------------
// Remove any degenerate points.
//------------------------------------------------------------------------------
template<int Dimension>
void
QuantPLC<Dimension>::removeDegeneracies() {
  const auto N = m_hashes.size();
  std::map<CoordHash, unsigned> seen;
  std::vector<CoordHash> new_hashes;
  std::vector<int> oldToNew(N, -1);
  new_hashes.reserve(N);
  unsigned newIndx = 0;
  for (auto i = 0; i < N; ++i) {
    const auto& h = m_hashes[i];
    auto it = seen.find(h);
    if (it != seen.end()) {
      oldToNew[i] = it->second;
    } else {
      seen.emplace(h, newIndx);
      oldToNew[i] = newIndx;
      new_hashes.push_back(h);
      newIndx++;
    }
  }
  m_hashes = std::move(new_hashes);
  m_points.clear();
  m_points.reserve(m_hashes.size());
  for (auto& h : m_hashes) {
    m_points.push_back(m_Q.unhash(h));
  }
  for (auto& f : facets) {
    for (auto& idx : f) {
      idx = oldToNew[idx];
    }
  }
  for (auto& hole : holes) {
    for (auto& f : hole) {
      for (auto& idx : f) {
        idx = oldToNew[idx];
      }
    }
  }
}
//------------------------------------------------------------------------------
// Remove any points not associated with a facet and order the facets.
//------------------------------------------------------------------------------
template<int Dimension>
void
QuantPLC<Dimension>::reduce() {
  // Collect only unique point indices used in the PLC facets
  std::set<int> indices;
  for (const auto& f : facets) {
    std::copy(f.begin(), f.end(), std::inserter(indices, indices.end()));
  }

  unsigned newIdx = 0;
  std::map<int, int> old2new;
  std::map<CoordHash, unsigned> hashToIndex;
  std::vector<CoordHash> newHashes;
  for(int oldIdx : indices) {
    auto hash = m_hashes[oldIdx];
    if (hashToIndex.find(hash) != hashToIndex.end()) {
      old2new[oldIdx] = hashToIndex[hash];
    } else {
      old2new[oldIdx] = newIdx;
      hashToIndex[hash] = newIdx++;
      newHashes.push_back(hash);
    }
  }

  auto N = indices.size();
  m_hashes = newHashes;
  m_points.clear();
  m_points.reserve(N);
  for (auto h : m_hashes) {
    m_points.push_back(m_Q.unhash(h));
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
  reduce();
}

// 2D
//---------------------
template<int Dimension>
template<int D>
std::enable_if_t<D == 2, void>
QuantPLC<Dimension>::makeConvex2D() {
  const unsigned n = m_points.size();
  m_convex_hull = true;

  // Unhash all points to integer coordinates and pair with original indices
  std::vector<std::pair<IntPoint, unsigned>> sortedPoints;
  sortedPoints.reserve(n);
  for (unsigned i = 0; i < n; ++i) {
    sortedPoints.push_back(std::make_pair(m_points[i], i));
  }

  // Sort by x-coordinate, then y-coordinate (using Point's operator<)
  std::sort(sortedPoints.begin(), sortedPoints.end(),
            [](const std::pair<IntPoint, unsigned>& a,
               const std::pair<IntPoint, unsigned>& b) {
              return a.first < b.first;
            });

  // Check if points are collinear
  bool collinear = true;
  if (n > 2) {
    for (unsigned i = 2; i < n && collinear; ++i) {
      collinear = (zcross_sign(sortedPoints[0].first,
                               sortedPoints[1].first,
                               sortedPoints[i].first) == 0);
    }
  }

  // If collinear, the hull is just a line segment from first to last
  if (collinear) {
    facets.resize(1, std::vector<int>(2));
    facets[0][0] = sortedPoints.front().second;
    facets[0][1] = sortedPoints.back().second;
    return;
  }

  // Non-collinear case: use monotone chain algorithm
  std::vector<int> result(2 * n);
  unsigned k = 0;

  // Build lower hull
  for (unsigned i = 0; i < n; ++i) {
    while (k >= 2 &&
           zcross_sign(sortedPoints[result[k - 2]].first,
                       sortedPoints[result[k - 1]].first,
                       sortedPoints[i].first) <= 0) {
      k--;
    }
    result[k++] = i;
  }

  // Build upper hull
  unsigned t = k + 1;
  for (int i = n - 2; i >= 0; --i) {
    while (k >= t &&
           zcross_sign(sortedPoints[result[k - 2]].first,
                       sortedPoints[result[k - 1]].first,
                       sortedPoints[i].first) <= 0) {
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
    facets.push_back(std::vector<int>(2));
    facets.back()[0] = sortedPoints[result[i]].second;
    facets.back()[1] = sortedPoints[result[i + 1]].second;
  }
}

// 3D
//---------------------
template<int Dimension>
template<int D>
std::enable_if_t<D == 3, void>
QuantPLC<Dimension>::makeConvex3D() {
  using RealPoint = Point<Dimension, double>;
  const unsigned n = m_points.size();
  m_convex_hull = true;

  std::vector<double> q_points;
  q_points.reserve(n*3);
  for (auto i = 0; i < n; ++i) {
    IntPoint ps = m_points[i];
    RealPoint rpoint = ps.template type_cast<double>();
    q_points.push_back(rpoint.x);
    q_points.push_back(rpoint.y);
    q_points.push_back(rpoint.z);
  }
  orgQhull::Qhull qh;
  qh.runQhull("", Dimension, n, q_points.data(), "Qt");

  facets.clear();
  facets.reserve(qh.facetList().size());
  // Extract and reformat the Qhull points
  for (auto f = qh.facetList().begin(); f != qh.facetList().end(); ++f) {
    auto vertices = f->vertices();
    facets.push_back(std::vector<int>(vertices.size()));
    unsigned fIndx = 0;
    for (auto v = vertices.begin(); v != vertices.end(); ++v, fIndx++) {
      auto pid = (*v).point().id();
      facets.back()[fIndx] = pid;
    }
  }
  // Remove any points that are not part of the convex hull
  reduce();
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
//------------------------------------------------------------------------------
template<int Dimension>
template<int D>
std::enable_if_t<D == 2, void>
QuantPLC<Dimension>::orderFacets2D() {
  auto orderEdgeLoop = [](std::vector<std::vector<int>>& edges) {
    if (edges.empty()) return;

    std::vector<std::vector<int>> ordered;
    ordered.reserve(edges.size());

    // Build map: start vertex -> edge index
    std::map<int, int> startMap;
    for (size_t i = 0; i < edges.size(); ++i) {
      startMap[edges[i][0]] = i;
    }

    // Follow the chain starting from first edge
    std::set<int> used;
    int current = 0;
    while (used.size() < edges.size()) {
      ordered.push_back(edges[current]);
      used.insert(current);

      int nextVertex = edges[current][1];
      if (startMap.count(nextVertex) && !used.count(startMap[nextVertex])) {
        current = startMap[nextVertex];
      } else {
        break;  // Chain broken
      }
    }

    edges = ordered;
  };

  orderEdgeLoop(facets);
  for (auto& hole : holes) {
    orderEdgeLoop(hole);
  }
}

//------------------------------------------------------------------------------
// Order 3D facets to have consistent winding
// Ensures all facet normals point outward from the centroid
// Works for arbitrary polygonal facets, not just triangles
//------------------------------------------------------------------------------
template<int Dimension>
template<int D>
std::enable_if_t<D == 3, void>
QuantPLC<Dimension>::orderFacets3D() {
  if (facets.empty()) return;
  using Wide = WideInt<3>;
  using WidePoint = Point3<Wide>;

  // Compute centroid
  const auto N = static_cast<Wide>(m_points.size());
  WidePoint centroid(0, 0, 0);
  for (const auto& p : m_points) {
    centroid = centroid + p.template type_cast<Wide>();
  }

  auto orientFacet = [this, &centroid, &N](std::vector<int>& facet) {
    if (facet.size() < 3) return;

    // Compute normal using first 3 vertices
    const IntPoint& v0 = m_points[facet[0]];
    const IntPoint& v1 = m_points[facet[1]];
    const IntPoint& v2 = m_points[facet[2]];
    auto normal = computeNormal(v0, v1, v2);

    // Vector from centroid to facet (use v0 as representative point)
    IntPoint toFace = normalizeByBitShift<IntType>(N*v0.template type_cast<Wide>() - centroid);
    Wide dot = qdot<3, IntType>(normal, toFace);

    // If dot < 0, normal points inward - reverse vertex order
    if (dot < 0) {
      std::reverse(facet.begin(), facet.end());
    }
  };

  for (auto& facet : facets) {
    orientFacet(facet);
  }

  for (auto& hole : holes) {
    for (auto& facet : hole) {
      orientFacet(facet);
      // For holes, we actually want inward normals, so reverse the result
      std::reverse(facet.begin(), facet.end());
    }
  }
}

//------------------------------------------------------------------------------
// Intersection tests
//------------------------------------------------------------------------------
template<int Dimension>
bool
QuantPLC<Dimension>::within(const RealPoint& point) const {
  IntPoint ip = m_Q.quantize(point);
  return within(ip);
}

template<int Dimension>
bool
QuantPLC<Dimension>::within(const IntPoint& point) const {
  POLY_ASSERT2(facets.size() > 0,
               "Must have facets before checking within");
  // First checks the bounds
  if (point < m_loBounds || point > m_hiBounds) {
    return false;
  }
  CoordHash phash = Hasher::hash(point);
  // Check if point is coincident with any vertices
  for (const auto& h : m_hashes) {
    if (h == phash) {
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
QuantPLC<Dimension>::within2D(const IntPoint& point) const {
  // Build the outer boundary polygon from ordered facets (edges)
  // Since facets are now ordered, simply collect first vertex of each edge
  std::vector<IntPoint> polygon;
  polygon.reserve(facets.size());
  for (const auto& f : facets) {
    polygon.push_back(m_points[f[0]]);
  }

  // Check if point is in or on the outer boundary
  if (polygon.size() < 3) return false;  // Degenerate case

  bool inOuter = pointInPolygon(point, polygon) || pointOnPolygon(point, polygon);
  if (!inOuter) return false;

  // Check if point is in any holes (if so, it's outside)
  for (const auto& hole : holes) {
    std::vector<IntPoint> holePoly;
    holePoly.reserve(hole.size());
    for (const auto& f : hole) {
      holePoly.push_back(m_points[f[0]]);
    }
    if (holePoly.size() >= 3 && pointInPolygon(point, holePoly)) {
      return false;
    }
  }

  return true;
}

template<int Dimension>
template<int D>
std::enable_if_t<D == 3, bool>
QuantPLC<Dimension>::within3D(const IntPoint& point) const {
  // Point is not on any boundary - use ray casting
  IntPoint rayEnd = point;
  rayEnd.x = m_hiBounds.x;

  // Check each face
  std::set<CoordHash> crossings;
  for (const auto& f : facets) {
    IntPoint hitPoint;
    int check = segmentFaceIntersection3D(point, rayEnd, f, m_points, hitPoint);
    if (check == 0) { // Coplanar and contained in a facet
      return true;
    } else if (check == 1) {
      crossings.insert(Hasher::hash(hitPoint));
    }
  }

  for (const auto& hole : holes) {
    for (const auto& f : hole) {
      IntPoint hitPoint;
      int check = segmentFaceIntersection3D(point, rayEnd, f, m_points, hitPoint);
      if (check == 0) {
        return true;
      } else if (check == 1) {
        crossings.insert(Hasher::hash(hitPoint));
      }
    }
  }
  // Odd number of crossings means inside
  return (crossings.size() % 2) == 1;
}

//------------------------------------------------------------------------------
// Instantiate versions we know we need.
//------------------------------------------------------------------------------
template class QuantPLC<2>;
template class QuantPLC<3>;
}
