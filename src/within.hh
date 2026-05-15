#ifndef POLYTOPE_WITHIN_HH
#define POLYTOPE_WITHIN_HH
//------------------------------------------------------------------------------
// within - determine whether a point lies inside a complex boundary
//
// Checks if a point is inside a PLC boundary AND outside holes, if present.
// By convention, if a point is on a hole boundary, it is classified as inside
// the full boundary.
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

// Determine if point lies inside a set of facets
// Functor definition first.
template<int Dimension, typename RealType> struct WithinBoundaryFunctor;

// 2D specialization.
template<typename RealType>
struct WithinBoundaryFunctor<2, RealType> {
  static int impl(const RealType* point,
		  const unsigned numVertices,
		  const RealType* vertices,
		  const std::vector<std::vector<int> >& facets) {
    const unsigned numFacets = facets.size();
    unsigned i = facets[0][0];

    // Check if point is on the boundary of the facet
    bool onBoundary = geometry::pointOnPolygon(point, numFacets, &vertices[2*i]);

    // Check if point is on the interior
    bool inBoundary = geometry::pointInPolygon(point, numFacets, &vertices[2*i]);

    if (onBoundary)                        return 2;
    else if(inBoundary and not onBoundary) return 1;
    else                                   return 0;
  }
};

// 3D specialization.
template<typename RealType>
struct WithinBoundaryFunctor<3, RealType> {
  static int impl(const RealType* point,
		  const unsigned numVertices,
		  const RealType* vertices,
		  const std::vector<std::vector<int> >& facets) {
    // Check if point is on the boundary of the polyhedron
    bool onBoundary = geometry::pointOnPolyhedron(point, facets, vertices);

    // Check if point is in the interior
    bool inBoundary = geometry::pointInPolyhedron(point, facets, vertices);

    if (onBoundary)                        return 2;
    else if(inBoundary and not onBoundary) return 1;
    else                                   return 0;
  }
};

// Functional interface.
template<int Dimension, typename RealType>
int withinBoundary(const RealType* point,
		   const unsigned numVertices,
		   const RealType* vertices,
		   const std::vector<std::vector<int> >& facets) {
  return WithinBoundaryFunctor<Dimension, RealType>::impl(point, numVertices, vertices, facets);
}

}

//------------------------------------------------------------------------------
// within
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
bool
within(const RealType* point,
       const unsigned numVertices,
       const RealType* vertices,
       const PLC<Dimension>& plc) {
  const int numHoles = plc.holes.size();
  int i = 0, insideType;

  // Check the outer boundary of the PLC.
  insideType = withinBoundary<Dimension, RealType>(point, numVertices, vertices, plc.facets);
  bool inBoundary = (insideType > 0) ? true : false;

  // Check each of the holes.
  bool inHole = false;
  while (i < numHoles and not inHole) {
    insideType = withinBoundary<Dimension, RealType>(point, numVertices, vertices, plc.holes[i]);
    inHole = (insideType == 1) ? true : false;
    ++i;
  }

  return (inBoundary and not inHole);
}

}

#endif
