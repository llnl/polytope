#ifndef __polytope_convexIntersection__
#define __polytope_convexIntersection__

#include <vector>
#include <utility>

#include "PLC.hh"

namespace { // anonymous

//------------------------------------------------------------------------------
// Compare a point to a line and determine if the point is above,
// below, or colinear with the line (assuming the line points (l1, l2) are
// specified in counter-clockwise manner to indicate the interior direction.
//------------------------------------------------------------------------------
template<typename RealType>
int compare(const RealType& l1x, const RealType& l1y,
            const RealType& l2x, const RealType& l2y,
            const RealType& px, const RealType& py) {
  const double ztest = (double(l2x - l1x)*double(py - l1y) -
                        double(l2y - l1y)*double(px - l1x));
  return -(ztest < 0.0 ? -1 :
           ztest > 0.0 ?  1 :
           0);
}

//------------------------------------------------------------------------------
// Compare a cloud of points to a line and determine if the points are above,
// below, or if the line passes through the points.  Assumes the line points 
// (l1, l2) are specified in counter-clockwise manner to indicate the interior 
// direction.
//------------------------------------------------------------------------------
template<typename RealType>
int compare(const RealType& l1x, const RealType& l1y, 
            const RealType& l2x, const RealType& l2y, 
            const std::vector<RealType>& points) {
  POLY_ASSERT(points.size() % 2 == 0);
  POLY_ASSERT(points.size() > 1);
  const unsigned n = points.size() / 2;
  const int result = compare(l1x, l1y, l2x, l2y, points[0], points[1]);
  unsigned i = 1;
  while (i < n and result == compare(l1x, l1y, l2x, l2y, 
                                     points[2*i], points[2*i + 1])) ++i;
  if (i == n) {
    return result;
  } else {
    return 0;
  }
}

//------------------------------------------------------------------------------
// Compare a point to a plane, and determine if the point is above,
// below, or coplanar with the plane.  The plane is specified by the 
// the point normal [(ox, oy, oz), (nx, ny, nz)], and the point by (px, py, pz).
//------------------------------------------------------------------------------
template<typename RealType>
int compare(const RealType& ox, const RealType& oy, const RealType& oz,
            const double& nx, const double& ny, const double& nz,
            const RealType& px, const RealType& py, const RealType& pz) {
  const double ztest = nx*double(px - ox) + ny*double(py - oy) + nz*double(pz - oz);
  return (ztest < 0.0 ? -1 :
          ztest > 0.0 ?  1 :
                         0);
}

//------------------------------------------------------------------------------
// Compare a cloud of points to a plane and determine if the points are above,
// below, or if the plane passes through the points.
//------------------------------------------------------------------------------
template<typename RealType>
int compare(const RealType& ox, const RealType& oy, const RealType& oz, 
            const RealType& nx, const RealType& ny, const RealType& nz, 
            const std::vector<RealType>& points) {
  POLY_ASSERT(points.size() % 3 == 0);
  POLY_ASSERT(points.size() > 1);
  const unsigned n = points.size() / 3;
  const int result = compare(ox, oy, oz, nx, ny, nz, points[0], points[1], points[2]);
  unsigned i = 1;
  while (i < n and result == compare(ox, oy, oz, nx, ny, nz,
                                     points[3*i], points[3*i + 1], points[3*i + 2])) ++i;
  if (i == n) {
    return result;
  } else {
    return 0;
  }
}

//------------------------------------------------------------------------------
// Compute the 3D normal given three points (a, b, c).
//------------------------------------------------------------------------------
template<typename RealType>
void computeNormal(const RealType& ax, const RealType& ay, const RealType& az,
                   const RealType& bx, const RealType& by, const RealType& bz,
                   const RealType& cx, const RealType& cy, const RealType& cz,
                   double& nx, double& ny, double& nz) {
  const double dx_ab = bx - ax, dy_ab = by - ay, dz_ab = bz - az;
  const double dx_ac = cx - ax, dy_ac = cy - ay, dz_ac = cz - az;
  nx = dy_ab*dz_ac - dz_ab*dy_ac;
  ny = dz_ab*dx_ac - dx_ab*dz_ac;
  nz = dx_ab*dy_ac - dy_ab*dx_ac;
}

//------------------------------------------------------------------------------
// Hash a pair points to represent an edge.
//------------------------------------------------------------------------------
std::pair<int, int>
hashEdge(const int i, const int j) {
  POLY_ASSERT(i != j);
  return (i < j ? std::make_pair(i, j) : std::make_pair(j, i));
}

}           // anonymous

namespace polytope {

//------------------------------------------------------------------------------
// Convex polygon intersection.
// We restrict this to ReducedPLC for expediency because the ReducedPLC has
// already computed the unique set of vertex coordinates.
// Vertex-to-edge, edge-to-edge, and vertex-to-vertex touching counts as intersection.
//------------------------------------------------------------------------------
template<typename RealType>
bool
convexIntersect(const ReducedPLC<2, RealType>& a, const ReducedPLC<2, RealType>& b) {
  const unsigned nva = a.points.size() / 2;
  const unsigned nvb = b.points.size() / 2;
  const unsigned nfa = a.facets.size();
  const unsigned nfb = b.facets.size();
  POLY_CONTRACT_VAR(nva);
  POLY_CONTRACT_VAR(nvb);
  unsigned i, j, k, m;

  // Any vertex containment is sufficient.
  for (i = 0; i < nva; ++i) {
    if (within<2, RealType>(&a.points[2*i], nvb, &b.points[0], b)) return true;
  }
  for (i = 0; i < nvb; ++i) {
    if (within<2, RealType>(&b.points[2*i], nva, &a.points[0], a)) return true;
  }

  // Otherwise, the polygons intersect only if an edge pair crosses.
  RealType intersectionPoint[2];
  for (i = 0; i < nfa; ++i) {
    j = a.facets[i][0];
    k = a.facets[i][1];
    POLY_ASSERT(j < nva);
    POLY_ASSERT(k < nva);
    for (m = 0; m < nfb; ++m) {
      const unsigned p = b.facets[m][0];
      const unsigned q = b.facets[m][1];
      POLY_ASSERT(p < nvb);
      POLY_ASSERT(q < nvb);
      if (geometry::segmentIntersection2D(&a.points[2*j], &a.points[2*k],
                                          &b.points[2*p], &b.points[2*q],
                                          intersectionPoint)) {
        return true;
      }
    }
  }

  return false;
}

//------------------------------------------------------------------------------
// Convex polygon intersection.
// This only works if the vertices are in counter-clockwise order.
//------------------------------------------------------------------------------
template<typename RealType>
bool
convexIntersectOrdered(const ReducedPLC<2, RealType>& a, const ReducedPLC<2, RealType>& b) {
  const unsigned nva = a.points.size() / 2;
  const unsigned nvb = b.points.size() / 2;
  const unsigned nfa = a.facets.size();
  const unsigned nfb = b.facets.size();
  POLY_CONTRACT_VAR(nva);
  POLY_CONTRACT_VAR(nvb);

  bool outside = false;
  unsigned i, j, ifacet;

  // Check if we can exclude b from a.
  {
    ifacet = 0;
    while (not outside and ifacet < nfa) {
      i = a.facets[ifacet][0];
      j = a.facets[ifacet][1];
      POLY_ASSERT(i < nva);
      POLY_ASSERT(j < nva);
      outside = (compare(a.points[2*i], a.points[2*i + 1],
                         a.points[2*j], a.points[2*j + 1],
                         b.points) == 1);
      ++ifacet;
    }
    if (outside) return false;
  }

  // Check if we can exclude a from b.
  {
    ifacet = 0;
    while (not outside and ifacet < nfb) {
      i = b.facets[ifacet][0];
      j = b.facets[ifacet][1];
      POLY_ASSERT(i < nvb);
      POLY_ASSERT(j < nvb);
      outside = (compare(b.points[2*i], b.points[2*i + 1],
                         b.points[2*j], b.points[2*j + 1],
                         a.points) == 1);
      ++ifacet;
    }
    if (outside) return false;
  }

  // We can't exclude anybody, so must intersect!
  return true;
}

//------------------------------------------------------------------------------
// Convex polyhedron intersection.
//------------------------------------------------------------------------------
template<typename RealType>
bool
convexIntersect(const ReducedPLC<3, RealType>& a, const ReducedPLC<3, RealType>& b) {
  const unsigned nva = a.points.size() / 3;
  const unsigned nvb = b.points.size() / 3;
  const unsigned nfa = a.facets.size();
  const unsigned nfb = b.facets.size();
  POLY_CONTRACT_VAR(nva);
  POLY_CONTRACT_VAR(nvb);

  bool outside = false;
  unsigned i, j, k, n, ifacet;
  double nx, ny, nz;

  // Check if we can exclude b from a.
  {
    ifacet = 0;
    while (not outside and ifacet < nfa) {
      i = a.facets[ifacet][0];
      j = a.facets[ifacet][1];
      k = a.facets[ifacet][2];
      POLY_ASSERT(i < nva);
      POLY_ASSERT(j < nva);
      POLY_ASSERT(k < nva);
      computeNormal(a.points[3*i], a.points[3*i + 1], a.points[3*i + 2],
                    a.points[3*j], a.points[3*j + 1], a.points[3*j + 2],
                    a.points[3*k], a.points[3*k + 1], a.points[3*k + 2],
                    nx, ny, nz);
      outside = (compare(a.points[3*i], a.points[3*i + 1], a.points[3*i + 2],
                         nx, ny, nz,
                         b.points) == 1);
      ++ifacet;
    }
    if (outside) return false;
  }

  // Check if we can exclude a from b.
  {
    ifacet = 0;
    while (not outside and ifacet < nfb) {
      i = b.facets[ifacet][0];
      j = b.facets[ifacet][1];
      k = b.facets[ifacet][2];
      POLY_ASSERT(i < nvb);
      POLY_ASSERT(j < nvb);
      POLY_ASSERT(k < nvb);
      computeNormal(b.points[3*i], b.points[3*i + 1], b.points[3*i + 2],
                    b.points[3*j], b.points[3*j + 1], b.points[3*j + 2],
                    b.points[3*k], b.points[3*k + 1], b.points[3*k + 2],
                    nx, ny, nz);
      outside = (compare(b.points[3*i], b.points[3*i + 1], b.points[3*i + 2],
                         nx, ny, nz,
                         a.points) == 1);
      ++ifacet;
    }
    if (outside) return false;
  }

  // Find the edges for each polyhedron.
  typedef std::pair<int, int> Edge;
  typedef std::set<Edge> EdgeSet;
  EdgeSet aEdges, bEdges;
  for (ifacet = 0; ifacet != nfa; ++ifacet) {
    n = a.facets[ifacet].size();
    for (i = 0; i != n; ++i) {
      j = (i + 1) % n;
      aEdges.insert(hashEdge(a.facets[ifacet][i], a.facets[ifacet][j]));
    }
  }
  for (ifacet = 0; ifacet != nfb; ++ifacet) {
    n = b.facets[ifacet].size();
    for (i = 0; i != n; ++i) {
      j = (i + 1) % n;
      bEdges.insert(hashEdge(b.facets[ifacet][i], b.facets[ifacet][j]));
    }
  }

  // Test against the cross products of the edges.
  int sidea, sideb;
  for (typename EdgeSet::const_iterator aItr = aEdges.begin();
       aItr != aEdges.end();
       ++aItr) {
    i = aItr->first;
    for (typename EdgeSet::const_iterator bItr = bEdges.begin();
         bItr != bEdges.end();
         ++bItr) {
      j = bItr->first;
      computeNormal(RealType(0), RealType(0), RealType(0),
                    a.points[3*i], a.points[3*i + 1], a.points[3*i + 2],
                    b.points[3*j], b.points[3*j + 1], b.points[3*j + 2],
                    nx, ny, nz);

      // Test all of a.
      sidea = compare(a.points[3*i], a.points[3*i + 1], a.points[3*i + 2],
                      nx, ny, nz,
                      a.points);
      if (sidea == 0) continue;
      sideb = compare(a.points[3*i], a.points[3*i + 1], a.points[3*i + 2],
                      nx, ny, nz,
                      b.points);
      if (sideb == 0) continue;
      if (sidea*sideb < 0) return false;
    }
  }

  // We can't exclude anybody, so must intersect!
  return true;
}

}

#endif
