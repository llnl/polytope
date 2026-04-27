#ifndef __polytope_convexIntersection__
#define __polytope_convexIntersection__

#include <vector>
#include <utility>

#include "PLC.hh"
#include "polytope_geometric_utilities.hh"

namespace { // anonymous

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
      outside = (geometry::aboveBelow(a.points[2*i], a.points[2*i + 1],
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
      outside = (geometry::aboveBelow(b.points[2*i], b.points[2*i + 1],
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
// Forward declaration for 3D ordered version.
//------------------------------------------------------------------------------
template<typename RealType>
bool convexIntersectOrdered(const ReducedPLC<3, RealType>& a, const ReducedPLC<3, RealType>& b);

//------------------------------------------------------------------------------
// Convex polyhedron intersection (ordering-independent).
//
// NOTE: This implementation requires 3D geometric primitives that are not yet
// implemented in polytope:
//   - pointInPolyhedron() / pointOnPolyhedron() (for vertex containment tests)
//   - segmentFaceIntersection3D() (for edge-face crossing tests)
//
// Until these are implemented, use convexIntersectOrdered() which requires
// properly oriented facets but uses the separating axis theorem.
//------------------------------------------------------------------------------
template<typename RealType>
bool
convexIntersect(const ReducedPLC<3, RealType>& a, const ReducedPLC<3, RealType>& b) {
  // For now, delegate to the ordered version which works if facets are properly oriented
  return convexIntersectOrdered(a, b);

  // TODO: Implement the ordering-independent version following this pattern:
  //
  // const unsigned nva = a.points.size() / 3;
  // const unsigned nvb = b.points.size() / 3;
  // const unsigned nfa = a.facets.size();
  // const unsigned nfb = b.facets.size();
  // POLY_CONTRACT_VAR(nva);
  // POLY_CONTRACT_VAR(nvb);
  // unsigned i, j, k, m;
  //
  // // Any vertex containment is sufficient.
  // for (i = 0; i < nva; ++i) {
  //   if (within<3, RealType>(&a.points[3*i], nvb, &b.points[0], b)) return true;
  // }
  // for (i = 0; i < nvb; ++i) {
  //   if (within<3, RealType>(&b.points[3*i], nva, &a.points[0], a)) return true;
  // }
  //
  // // Otherwise, the polyhedra intersect only if an edge crosses a face.
  // // For each edge of a, test against each face of b.
  // for (i = 0; i < nfa; ++i) {
  //   const unsigned nEdges = a.facets[i].size();
  //   for (j = 0; j < nEdges; ++j) {
  //     k = (j + 1) % nEdges;
  //     const unsigned e0 = a.facets[i][j];
  //     const unsigned e1 = a.facets[i][k];
  //     POLY_ASSERT(e0 < nva and e1 < nva);
  //     for (m = 0; m < nfb; ++m) {
  //       if (geometry::segmentFaceIntersection3D(&a.points[3*e0], &a.points[3*e1],
  //                                               b.facets[m], &b.points[0])) {
  //         return true;
  //       }
  //     }
  //   }
  // }
  //
  // // For each edge of b, test against each face of a.
  // for (i = 0; i < nfb; ++i) {
  //   const unsigned nEdges = b.facets[i].size();
  //   for (j = 0; j < nEdges; ++j) {
  //     k = (j + 1) % nEdges;
  //     const unsigned e0 = b.facets[i][j];
  //     const unsigned e1 = b.facets[i][k];
  //     POLY_ASSERT(e0 < nvb and e1 < nvb);
  //     for (m = 0; m < nfa; ++m) {
  //       if (geometry::segmentFaceIntersection3D(&b.points[3*e0], &b.points[3*e1],
  //                                               a.facets[m], &a.points[0])) {
  //         return true;
  //       }
  //     }
  //   }
  // }
  //
  // return false;
}

//------------------------------------------------------------------------------
// Convex polyhedron intersection (ordered version).
// This only works if the facet vertices are properly oriented (outward-facing
// normals) using the right-hand rule. Uses the separating axis theorem.
//------------------------------------------------------------------------------
template<typename RealType>
bool
convexIntersectOrdered(const ReducedPLC<3, RealType>& a, const ReducedPLC<3, RealType>& b) {
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
      geometry::computeNormal(a.points[3*i], a.points[3*i + 1], a.points[3*i + 2],
                              a.points[3*j], a.points[3*j + 1], a.points[3*j + 2],
                              a.points[3*k], a.points[3*k + 1], a.points[3*k + 2],
                              nx, ny, nz);
      outside = (geometry::aboveBelow(a.points[3*i], a.points[3*i + 1], a.points[3*i + 2],
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
      geometry::computeNormal(b.points[3*i], b.points[3*i + 1], b.points[3*i + 2],
                              b.points[3*j], b.points[3*j + 1], b.points[3*j + 2],
                              b.points[3*k], b.points[3*k + 1], b.points[3*k + 2],
                              nx, ny, nz);
      outside = (geometry::aboveBelow(b.points[3*i], b.points[3*i + 1], b.points[3*i + 2],
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
      geometry::computeNormal(RealType(0), RealType(0), RealType(0),
                              a.points[3*i], a.points[3*i + 1], a.points[3*i + 2],
                              b.points[3*j], b.points[3*j + 1], b.points[3*j + 2],
                              nx, ny, nz);

      // Test all of a.
      sidea = geometry::aboveBelow(a.points[3*i], a.points[3*i + 1], a.points[3*i + 2],
                                   nx, ny, nz,
                                   a.points);
      if (sidea == 0) continue;
      sideb = geometry::aboveBelow(a.points[3*i], a.points[3*i + 1], a.points[3*i + 2],
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
