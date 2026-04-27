#ifndef POLYTOPE_INTERSECT_HH
#define POLYTOPE_INTERSECT_HH
//------------------------------------------------------------------------------
// intersect - compute the number of intersections of a line segment and a 
// PLC boundary
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
template<int Dimension, typename RealType> struct IntersectFacetsFunctor;

// 2D specialization.
template<typename RealType>
struct IntersectFacetsFunctor<2, RealType> {
  static unsigned impl(const RealType* point1,
                       const RealType* point2,
                       const unsigned numVertices,
                       const RealType* vertices,
                       const std::vector<std::vector<int> >& facets,
                       std::vector<RealType>& result) {
    unsigned i, j, numIntersections=0;
    RealType intersectionPoint[2];
    bool intersects, addPoint;
    const unsigned numFacets = facets.size();
    for (unsigned ifacet = 0; ifacet != numFacets; ++ifacet) {
      POLY_ASSERT(facets[ifacet].size() == 2);
      i = facets[ifacet][0];
      j = facets[ifacet][1];
      // POLY_ASSERT(i >= 0 and i < numVertices);
      // POLY_ASSERT(j >= 0 and j < numVertices);
      intersects = geometry::segmentIntersection2D(point1, point2,
                                                   &vertices[2*i], &vertices[2*j],
                                                   intersectionPoint);
      addPoint = false;
      if(intersects) {
        if (numIntersections == 0) addPoint = true;
        else {
          if (intersectionPoint[0] != result.back() - 1 and
              intersectionPoint[1] != result.back() ) addPoint = true;
          else addPoint = false;
        }
      }
      
      if (addPoint) {
        ++numIntersections;
        result.push_back(intersectionPoint[0]);
        result.push_back(intersectionPoint[1]);
      }
    }
    return numIntersections;
  }
};

// Functional interface.
template<int Dimension, typename RealType> 
unsigned intersectFacets(const RealType* point1,
                         const RealType* point2,
                         const unsigned numVertices,
                         const RealType* vertices,
                         const std::vector<std::vector<int> >& facets,
                         std::vector<RealType>& result) {
  return IntersectFacetsFunctor<2, RealType>::impl(point1, point2, numVertices, vertices, facets, result);
}

}


//------------------------------------------------------------------------------
// intersect
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
unsigned
intersect(const RealType* point1,
          const RealType* point2,
          const unsigned numVertices,
          const RealType* vertices,
          const PLC<Dimension, RealType>& plc,
          std::vector<RealType>& result) {

  // Check the outer boundary of the PLC.
  unsigned numIntersections = intersectFacets<Dimension, RealType>(point1, point2, numVertices, vertices, plc.facets, result);

  // Check each of the holes.
  unsigned numHoles = plc.holes.size();
  std::vector<unsigned> holeIntersections(numHoles);
  for (unsigned ihole = 0; ihole != plc.holes.size(); ++ihole) {
    holeIntersections[ihole] = intersectFacets<Dimension, RealType>(point1, point2, numVertices, vertices, plc.holes[ihole], result);
    numIntersections += holeIntersections[ihole];
  }
  return numIntersections;
}

}

#endif
