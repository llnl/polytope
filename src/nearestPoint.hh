#ifndef POLYTOPE_NEAREST_POINT_HH
#define POLYTOPE_NEAREST_POINT_HH
//------------------------------------------------------------------------------
// nearestPoint - compute the closest point on a PLC to a given point.
//
// The point can be on the interior or exterior of the PLC.  Checking of holes
// is also included.
// The point is returned as the last argument, and the distance as the result.
//------------------------------------------------------------------------------
#include <vector>
#include <limits>

#include "PLC.hh"
#include "polytope_internal.hh"
#include "polytope_geometric_utilities.hh"

namespace polytope {

//------------------------------------------------------------------------------
// Local utility methods.
//------------------------------------------------------------------------------
namespace {

// Find the closest point on a given set of facets.
// Functor definition first.
template<int Dimension, typename RealType> struct ClosestPointOnFacetsFunctor;

// 2D specialization.
template<typename RealType>
struct ClosestPointOnFacetsFunctor<2, RealType> {
  static RealType impl(const RealType* point,
                       const unsigned numVertices,
                       const RealType* vertices,
                       const std::vector<std::vector<int> >& facets,
                       RealType* result) {
    unsigned i, j;
    RealType dist, minDist = std::numeric_limits<RealType>::max();
    RealType candidate[2];
    (void) numVertices;
    const unsigned numFacets = facets.size();
    for (unsigned ifacet = 0; ifacet != numFacets; ++ifacet) {
      POLY_ASSERT(facets[ifacet].size() == 2);
      i = facets[ifacet][0];
      j = facets[ifacet][1];
      POLY_ASSERT(i < numVertices);
      POLY_ASSERT(j < numVertices);
      geometry::closestPointOnSegment2D(point, &vertices[2*i], &vertices[2*j], candidate);
      dist = geometry::distance<2, RealType>(point, candidate);
      if (dist < minDist) {
        minDist = dist;
        result[0] = candidate[0];
        result[1] = candidate[1];
      }
    }
    return minDist;
  }
};

// Functional interface.
template<int Dimension, typename RealType> 
RealType closestPointOnFacets(const RealType* point,
                              const unsigned numVertices,
                              const RealType* vertices,
                              const std::vector<std::vector<int> >& facets,
                              RealType* result) {
  return ClosestPointOnFacetsFunctor<2, RealType>::impl(point, numVertices, vertices, facets, result);
}

}

//------------------------------------------------------------------------------
// nearestPoint
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
RealType
nearestPoint(const RealType* point,
             const unsigned numVertices,
             const RealType* vertices,
             const PLC<Dimension>& plc,
             RealType* result) {

  // Check the outer boundary of the PLC.
  RealType minDist = closestPointOnFacets<Dimension, RealType>(point, numVertices, vertices, plc.facets, result);

  // Check each of the holes.
  RealType dist, candidate[2];
  for (unsigned ihole = 0; ihole != plc.holes.size(); ++ihole) {
    dist = closestPointOnFacets<Dimension, RealType>(point, numVertices, vertices, plc.holes[ihole], candidate);
    if (dist < minDist) {
      minDist = dist;
      for (unsigned i = 0; i != Dimension; ++i) result[i] = candidate[i];
    }
  }
  return minDist;
}

}

#endif
