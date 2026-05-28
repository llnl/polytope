//------------------------------------------------------------------------
// BoostTessellator
//------------------------------------------------------------------------
#include "BoostTessellator.hh"

#include <iostream>

// Handy Boost stuff
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "polytope.hh" // Pulls in POLY_ASSERT

// The Voronoi tools in Boost.Polygon
#include <boost/polygon/voronoi.hpp>

namespace polytope {

//------------------------------------------------------------------------------
// Compute the QuantizedTessellation
//------------------------------------------------------------------------------
template<typename RealType>
void
BoostTessellator<RealType>::
tessellateQuantized(QuantizedTessellation& result) const {

  // Get the quantized generator points as reals for Boost
  std::vector<Point2<RealType>> generators = result.getRealPoints();
  const size_t numGenerators = generators.size();

  // Invoke the Boost.Voronoi diagram constructor
  VD voronoi;
  construct_voronoi(generators.begin(), generators.end(), &voronoi);

  // Build the tessellation data structures
  // In 2D: nodes are Voronoi vertices, faces are edges, cells are Voronoi cells
  result.m_nodes.reserve(voronoi.num_vertices());
  result.m_faces.reserve(voronoi.num_edges());
  result.m_cells.resize(numGenerators);

  // Map Voronoi vertices to our node indices
  // Type aliases
  using IntType = typename QuantTessellation<2>::IntType;
  using IntPoint = typename QuantTessellation<2>::IntPoint;
  std::map<const typename VD::vertex_type*, int> vertex2node;

  // Process each Voronoi cell
  for (typename VD::const_cell_iterator cellItr = voronoi.cells().begin();
       cellItr != voronoi.cells().end();
       ++cellItr) {

    if (!cellItr->contains_point()) continue;

    const int cellIndex = cellItr->source_index();
    if (cellIndex >= numGenerators) continue;

    // Walk edges CCW around this cell
    const typename VD::edge_type* firstEdge = cellItr->incident_edge();
    const typename VD::edge_type* edge = firstEdge;

    do {
      const typename VD::vertex_type* v0 = edge->vertex0();
      const typename VD::vertex_type* v1 = edge->vertex1();

      // Skip infinite edges (should not happen with proper bounding)
      POLY_ASSERT2(v0 && v1, "Infinite edge encountered for cell " << cellIndex);

      // Add vertices to node list if not already present
      // Boost returns floating-point coordinates, but they're in quantized space
      // Round back to integers
      if (vertex2node.find(v0) == vertex2node.end()) {
        vertex2node[v0] = result.m_nodes.size();
        IntType x0 = static_cast<IntType>(std::round(v0->x()));
        IntType y0 = static_cast<IntType>(std::round(v0->y()));
        result.m_nodes.push_back(IntPoint(x0, y0));
      }
      if (vertex2node.find(v1) == vertex2node.end()) {
        vertex2node[v1] = result.m_nodes.size();
        IntType x1 = static_cast<IntType>(std::round(v1->x()));
        IntType y1 = static_cast<IntType>(std::round(v1->y()));
        result.m_nodes.push_back(IntPoint(x1, y1));
      }

      int n0 = vertex2node[v0];
      int n1 = vertex2node[v1];

      // Skip degenerate edges
      if (n0 != n1) {
        // Create face (edge in 2D) as a vector of 2 node indices
        int faceIndex = result.m_faces.size();
        result.m_faces.push_back({n0, n1});
        result.m_cells[cellIndex].push_back(faceIndex);
      }

      edge = edge->next();
    } while (edge != firstEdge);
  }
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class BoostTessellator<double>;

} //end polytope namespace
