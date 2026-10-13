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
VoronoiAssembler<2>
BoostTessellator::
tessellateQuantizedImpl(QuantizedTessellation& result) const {
  // Type aliases
  using VD = boost::polygon::voronoi_diagram<RealType>;
  const auto& Q = Quantizer<2>::instance();
  // First ensure that the encoding method has not changed
  POLY_CHECK2(Q.keyEncoding() == result.keyEncoding(),
              "Key encoding method changed during tessellation");
  // Get the quantized generators
  std::vector<QuantizedPoint<2>> generators = result.getQuantizedPoints();

  VD voronoi;

  // Invoke the Boost.Voronoi diagram constructor
  typedef boost::polygon::detail::voronoi_ctype_traits<QuantizedCoordinate<2>> MyTraits;
  boost::polygon::voronoi_builder<QuantizedCoordinate<2>, MyTraits> builder;
  for (const auto& p : generators) {
    builder.insert_point(p.x, p.y);
  }
  builder.construct(&voronoi);

  VoronoiAssembler<2> assembler(result);

  // Process each Voronoi edge
  for (const auto& edge : voronoi.edges()) {
    if (&edge > edge.twin()) {
      continue;
    }
    const auto* cell0 = edge.cell();
    const auto* cell1 = edge.twin()->cell();

    const auto& gindx0 = cell0->source_index();
    const auto& gindx1 = cell1->source_index();
    const auto* v0 = edge.vertex0();
    const auto* v1 = edge.vertex1();
    if (v0 && v1) {
      assembler.addFiniteEdge(gindx0, gindx1,
                              Point2<double>(v0->x(), v0->y()),
                              Point2<double>(v1->x(), v1->y()));
    } else if (v0) {
      assembler.addStartRay(gindx0, gindx1,
                            Point2<double>(v0->x(), v0->y()));
    } else if (v1) {
      assembler.addEndRay(gindx0, gindx1,
                          Point2<double>(v1->x(), v1->y()));
    } else {
      assembler.addInfLines(gindx0, gindx1);
    }
  }
  return assembler;
}

} //end polytope namespace
