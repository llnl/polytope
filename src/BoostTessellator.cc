//------------------------------------------------------------------------
// BoostTessellator
//------------------------------------------------------------------------
#include "BoostTessellator.hh"
#include "EdgeUtils.hh"

#include <iostream>
#include <utility>

// Handy Boost stuff
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "polytope_internal.hh" // Pulls in POLY_ASSERT
#include "RegisterBoostPolygonTypes.hh"

// The Voronoi tools in Boost.Polygon
#include <boost/polygon/voronoi.hpp>

namespace polytope {

//------------------------------------------------------------------------------
// Compute the QuantizedTessellation
//------------------------------------------------------------------------------
BoostTessellator::PrimitiveCells
BoostTessellator::
tessellateQuantizedImpl(const QuantizedTessellation& result) const {
  // Type aliases
  using VD = boost::polygon::voronoi_diagram<RealType>;
  const auto& Q = Quantizer<2>::instance();
  // First ensure that the encoding method has not changed
  POLY_CHECK2(Q.keyEncoding() == result.keyEncoding(),
              "Key encoding method changed during tessellation");
  // Get the quantized generators
  std::vector<QuantizedPoint<2>> generators = result.getQuantizedPoints();
  const size_t numGenerators = generators.size();

  VD voronoi;

  // Invoke the Boost.Voronoi diagram constructor
  typedef boost::polygon::detail::voronoi_ctype_traits<QuantizedCoordinate<2>> MyTraits;
  boost::polygon::voronoi_builder<QuantizedCoordinate<2>, MyTraits> builder;
  for (const auto& p : generators) {
    builder.insert_point(p.x, p.y);
  }
  builder.construct(&voronoi);

  PrimitiveCells cellPrimitives(numGenerators);

  // Process each Voronoi cell
  for (typename VD::const_cell_iterator cellItr = voronoi.cells().begin();
       cellItr != voronoi.cells().end();
       ++cellItr) {

    if (!cellItr->contains_point()) continue;

    const int cellIndex = cellItr->source_index();
    if (cellIndex >= int(numGenerators)) continue;

    // Walk edges CCW around this cell
    const typename VD::edge_type* firstEdge = cellItr->incident_edge();
    const typename VD::edge_type* edge = firstEdge;
    std::vector<VoronoiPrimitive<2>> vps;
    do {
      const VD::edge_type* nextEdge = edge->next();
      const typename VD::vertex_type* v0 = edge->vertex0();
      const typename VD::vertex_type* v1 = edge->vertex1();

      // An edge is considered infinite if Boost provides a null pointer
      // gen0 should always be the current cell's generator
      auto gindx1 = edge->cell()->source_index();
      auto gindx2 = edge->twin()->cell()->source_index();
      VoronoiPrimitive<2> vp(gindx1, gindx2);
      if (v0) {
        vp.setV0(Point2<double>(v0->x(), v0->y()));
      }
      if (v1) {
        vp.setV1(Point2<double>(v1->x(), v1->y()));
      }
      vps.push_back(vp);
      edge = nextEdge;
    } while (edge != firstEdge);
    cellPrimitives[cellIndex] = std::move(vps);
  }
  return cellPrimitives;
}

} //end polytope namespace
