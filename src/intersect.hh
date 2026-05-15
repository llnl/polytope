#ifndef POLYTOPE_INTERSECT_HH
#define POLYTOPE_INTERSECT_HH
//------------------------------------------------------------------------------
// intersect - compute the number of intersections of a line segment and a
// PLC boundary
//------------------------------------------------------------------------------
#include <vector>
#include <set>
#include <array>
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
  static std::set<std::array<RealType, 2>> impl(const RealType* point1,
                                                const RealType* point2,
                                                const unsigned numVertices,
                                                const RealType* vertices,
                                                const std::vector<std::vector<int> >& facets) {
    std::set<std::array<RealType, 2>> result;
    unsigned i, j;
    RealType intersectionPoint[2];
    const unsigned numFacets = facets.size();
    for (unsigned ifacet = 0; ifacet < numFacets; ++ifacet) {
      POLY_ASSERT(facets[ifacet].size() == 2);
      i = facets[ifacet][0];
      j = facets[ifacet][1];
      // POLY_ASSERT(i >= 0 and i < numVertices);
      // POLY_ASSERT(j >= 0 and j < numVertices);
      bool intersects = geometry::segmentIntersection2D(point1, point2,
                                                        &vertices[2*i], &vertices[2*j],
                                                        intersectionPoint);
      if (intersects) {
        // std::set automatically handles duplicates
        result.insert({intersectionPoint[0], intersectionPoint[1]});
      }
    }
    return result;
  }
};

// 3D specialization.
template<typename RealType>
struct IntersectFacetsFunctor<3, RealType> {
  static std::set<std::array<RealType, 3>> impl(const RealType* point1,
                                                const RealType* point2,
                                                const unsigned numVertices,
                                                const RealType* vertices,
                                                const std::vector<std::vector<int> >& facets) {
    std::set<std::array<RealType, 3>> result;
    const unsigned numFacets = facets.size();
    const RealType tol = 1.0e-10;
    RealType intersectionPoint[3];

    for (unsigned ifacet = 0; ifacet < numFacets; ++ifacet) {
      // Use the compute function to get the intersection point
      bool intersects = geometry::segmentFaceIntersection3D(
        point1, point2, facets[ifacet], vertices, intersectionPoint, tol);

      if (intersects) {
        // std::set automatically handles duplicates
        result.insert({intersectionPoint[0], intersectionPoint[1], intersectionPoint[2]});
      }
    }
    return result;
  }
};

// Functional interface.
template<int Dimension, typename RealType>
std::set<std::array<RealType, Dimension>> intersectFacets(const RealType* point1,
                                                          const RealType* point2,
                                                          const unsigned numVertices,
                                                          const RealType* vertices,
                                                          const std::vector<std::vector<int> >& facets) {
  return IntersectFacetsFunctor<Dimension, RealType>::impl(point1, point2, numVertices, vertices, facets);
}

}


//------------------------------------------------------------------------------
// Intersect
// Gather unique intersection points from segment between point1 and point2.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
std::set<std::array<RealType, Dimension>>
intersect(const RealType* point1,
          const RealType* point2,
          const unsigned numVertices,
          const RealType* vertices,
          const PLC<Dimension>& plc) {

  // Check the outer boundary of the PLC.
  auto result = intersectFacets<Dimension, RealType>(point1, point2, numVertices, vertices, plc.facets);

  // Check each of the holes.
  for (unsigned ihole = 0; ihole < plc.holes.size(); ++ihole) {
    auto holeIntersections = intersectFacets<Dimension, RealType>(point1, point2, numVertices, vertices, plc.holes[ihole]);
    result.insert(holeIntersections.begin(), holeIntersections.end());
  }

  return result;
}

}

#endif
