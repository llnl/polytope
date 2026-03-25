//------------------------------------------------------------------------
// BoostTessellator
//------------------------------------------------------------------------
#include "BoostTessellator.hh"

#include <iostream>
#include <algorithm>
#include <set>
#include <list>
#include <map>
#include <limits>
#include "float.h"

// Handy Boost stuff
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "polytope.hh" // Pulls in POLY_ASSERT

#include "within.hh"
#include "nearestPoint.hh"
#include "intersect.hh"
#include "polytope_plc_canned_geometries.hh"
#include "PLC_Boost_2d.hh"
#include "Segment.hh"

#include "RegisterBoostPolygonTypes.hh"

// Fast predicate for determining colinearity of points.
extern double orient2d(double* pa, double* pb, double* pc);

// The Voronoi tools in Boost.Polygon
#include <boost/polygon/voronoi.hpp>

namespace polytope {

using namespace std;
using std::min;
using std::max;
using std::abs;

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
template<typename RealType>
BoostTessellator<RealType>::
BoostTessellator():
  Tessellator<2, RealType>() {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename RealType>
BoostTessellator<RealType>::
~BoostTessellator() {
}

//------------------------------------------------------------------------------
// Compute the QuantizedTessellation
//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
tessellateQuantized(QuantizedTessellation& result) const {

  // Sort the input points by the first element of the generator-index pair.
  // The second element provides the pre-sort generator index. Shame on you,
  // Boost, for making us do this.
  const int numGenerators = result.generators.size();
  vector<IntPoint> generators(result.generators.begin(), result.generators.end());
  std::copy(result.guardGenerators.begin(), result.guardGenerators.end(), std::back_inserter(generators));
  // sort(generators.begin(), generators.end());
  
  // // Build ourselves the segments representing our bounding box.
  // typedef Segment<2> IntSegment;
  // vector<IntSegment> bounds(4);
  // const CoordHash coordMin = QuantizedTessellation2d<CoordHash, RealType>::coordMin;
  // const CoordHash coordMax = QuantizedTessellation2d<CoordHash, RealType>::coordMax;
  // bounds[0] = IntSegment(coordMin, coordMin, coordMax, coordMin);
  // bounds[1] = IntSegment(coordMax, coordMin, coordMax, coordMax);
  // bounds[2] = IntSegment(coordMax, coordMax, coordMin, coordMax);
  // bounds[3] = IntSegment(coordMin, coordMax, coordMin, coordMin);
  // // sort(bounds.begin(), bounds.end());  // Not sure if this is necessary

  // Invoke the Boost.Voronoi diagram constructor
  VD voronoi;
  construct_voronoi(generators.begin(), generators.end(),
                    // bounds.begin(), bounds.end(),
                    &voronoi);  
  // POLY_ASSERT(voronoi.num_cells() == numGenerators + result.guardGenerators.size());

  // Read out the Voronoi topology to our intermediate QuantizedTessellation format.
  // Iterate over the edges.  Boost has organized them CCW around each generator.
  // Note here we just extract information about the point generators, and build straight
  // edges where those interface with our segment boundaries.
  result.cellEdges = vector<vector<int> >(numGenerators);
  result.edges.reserve(voronoi.num_edges());
  result.nodes.reserve(voronoi.num_vertices());
  map<IntPoint, int> node2id;
  map<std::pair<int, int>, int> edge2id;
  for (typename VD::const_cell_iterator cellItr = voronoi.cells().begin(); 
       cellItr != voronoi.cells().end(); 
       ++cellItr) {
    int cellIndex = generators[cellItr->source_index()].index;
    if (cellItr->contains_point() and cellIndex < numGenerators) {
    // if (cellItr->contains_point() and cellItr->source_index() < numGenerators) {
      // POLY_ASSERT2(sortedIndex == cellItr->source_index(),
      //              sortedIndex << " != " << cellItr->source_index());
      // POLY_ASSERT(sortedIndex <  numGenerators + 4);
      // POLY_ASSERT(cellItr->source_index() <  numGenerators);
      POLY_ASSERT2(cellIndex   <  numGenerators, cellIndex << " " << cellItr->source_index() << " " << numGenerators);

      // Start the chain walking the edges of this cell.
      const typename VD::edge_type* edge = cellItr->incident_edge();
      do {
        edge = edge->next();

        // The two vertex pointers for this edge
        // NOTE: If edge is infinite, one of these pointers is null.  This should never happen
        // since we added bounding segments.
        const typename VD::vertex_type* v0 = edge->vertex0();
        const typename VD::vertex_type* v1 = edge->vertex1();
        if (not (v0 and v1)) {
          cerr << "Bad news at generator " << generators[cellIndex] << " " << result.coordMin << " " << result.coordMax << endl;
          if (v0) {
            cerr << "v0 : " << v0->x() << " " << v0->y() << endl;
          } else {
            cerr << "v1 : " << v1->x() << " " << v1->y() << endl;
          }
        }
        POLY_ASSERT(v0 and v1);

        // Finite edge.  Add the edge to the cell, and any new nodes.
        // Insert vertex 0.
        const IntPoint p0(v0->x(), v0->y());
        int old_size = node2id.size();
        const int j0 = internal::addKeyToMap(p0, node2id);
        if (j0 == old_size) {
          POLY_ASSERT(j0 == result.nodes.size());
          result.nodes.push_back(p0);
        }

        // Insert vertex 1.
        const IntPoint p1(v1->x(), v1->y());
        old_size = node2id.size();
        const int j1 = internal::addKeyToMap(p1, node2id);
        if (j1 == old_size) {
          POLY_ASSERT(j1 == result.nodes.size());
          result.nodes.push_back(p1);
        }
        
        // Now insert the edge.  Since we use oriented single edges between cells, a bit different than Boost.Polygon.
        // We have to screen out zero-length edges apparently.
        if (j0 != j1) {
          const pair<int, int> edge = internal::hashEdge(j0, j1);
          POLY_ASSERT((edge.first == j0 and edge.second == j1) or
                      (edge.first == j1 and edge.second == j0));
          old_size = edge2id.size();
          const int e1 = internal::addKeyToMap(edge, edge2id);
          if (e1 == old_size) {
            POLY_ASSERT(e1 == result.edges.size());
            result.edges.push_back(edge);
          }
          if (edge.first == j0) {
            result.cellEdges[cellIndex].push_back(e1);
          } else {
            result.cellEdges[cellIndex].push_back(~e1);
          }
        }
      } while (edge != cellItr->incident_edge());
    }
  }
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class BoostTessellator<double>;

} //end polytope namespace
