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
QuantPLC<Dimension>::
QuantPLC(const PLC<Dimension> plc,
         const Quant& Q,
         const std::vector<RealType>& allpoints,
         bool doReduce) :
  PLC<Dimension>(plc),
  m_Q(Q) {

  m_loBounds.one();
  m_hiBounds.zero();

  // Extract the unrolled coordinates
  std::vector<RealPoint> rpoints = extractCoords<Dimension, RealType>(allpoints);

  // Apply reduction operation
  if (doReduce && facets.size() > 0) {
    reduce(rpoints);
    m_reduced = true;
  } else {
    auto N = rpoints.size();
    m_hashes.reserve(N);
    m_points.reserve(N);
    for(const auto& rp : rpoints) {
      auto ip = m_Q.quantize(rp);
      m_loBounds = m_loBounds.minElements(ip);
      m_hiBounds = m_hiBounds.maxElements(ip);
      m_hashes.push_back(m_Q.hash(ip));
      m_points.push_back(ip);
    }
  }
  POLY_ASSERT2(m_loBounds < m_hiBounds,
               "Provided coplanar or collinear or degenerate points to the QuantPLC");
}

template<int Dimension>
void
QuantPLC<Dimension>::
reduce(const std::vector<RealPoint>& rpoints) {
  auto N = rpoints.size();
  m_points.clear();
  m_points.reserve(N);
  m_hashes.clear();
  m_hashes.reserve(N);
  m_loBounds.one();
  m_hiBounds.zero();
  // Collect only unique point indices used in the PLC facets
  std::set<int> indices;
  for (auto i = 0; i < facets.size(); ++i) {
    std::copy(facets[i].begin(), facets[i].end(),
              std::inserter(indices, indices.end()));
  }

  unsigned newIdx = 0;
  std::map<int, int> old2new;
  std::map<CoordHash, unsigned> hashToIndex;
  for(int oldIdx : indices) {
    auto ip = m_Q.quantize(rpoints[oldIdx]);
    m_loBounds = m_loBounds.minElements(ip);
    m_hiBounds = m_hiBounds.maxElements(ip);
    auto hash = m_Q.hash(ip);
    if (hashToIndex.find(hash) != hashToIndex.end()) {
      old2new[oldIdx] = hashToIndex[hash];
    } else {
      old2new[oldIdx] = newIdx;
      hashToIndex[hash] = newIdx;
      m_hashes.push_back(hash);
      m_points.push_back(ip);
      ++newIdx;
    }
  }

  // Copy facets with remapped indices
  for (size_t i = 0; i < facets.size(); ++i) {
    for (size_t j = 0; j < facets[i].size(); ++j) {
      facets[i][j] = old2new[facets[i][j]];
    }
  }
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
  // First checks the bounds
  if (point < m_loBounds || m_hiBounds < point) {
    return false;
  }
  CoordHash phash = Hasher::hash(point);
  // Check if point is coincident with any vertices
  for (const auto& hash : m_hashes) {
    if (hash == phash) {
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
  for (const auto& f : facets) {
    if (!pointInPolygon(point, f, m_points)) {
      return false;
    }
  }
  for (const auto& hole : holes) {
    for (const auto& f : hole) {
      if (pointInPolygon(point, f, m_points)) {
        return false;
      }
    }
  }
  return true;
}

template<int Dimension>
template<int D>
std::enable_if_t<D == 3, bool>
QuantPLC<Dimension>::within3D(const IntPoint& point) const {
  IntPoint rayEnd = point;
  rayEnd.x = m_hiBounds.x;
  // Check each face
  std::set<CoordHash> crossings;
  for (const auto& f : facets) {
    IntPoint hitPoint;
    if (segmentFaceIntersection3D(point, rayEnd, f, m_points, hitPoint)) {
      crossings.insert(Hasher::hash(hitPoint));
    }
  }

  for (const auto& hole : holes) {
    for (const auto& f : hole) {
      IntPoint hitPoint;
      if (segmentFaceIntersection3D(point, rayEnd, f, m_points, hitPoint)) {
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
