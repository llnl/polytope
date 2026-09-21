//------------------------------------------------------------------------
// TriangleTessellator
//
// Polytope wrapper for the native 2D Voronoi tessellator in Triangle.Polygon
// v1.52 or greater
//------------------------------------------------------------------------
#ifndef __Polytope_TriangleTessellator__
#define __Polytope_TriangleTessellator__

#include <vector>
#include <cmath>
#include <limits>

#include "SerialTessellator.hh"
#include "QuantTessellation.hh"
#include "Point.hh"
#include "VoronoiAssembler.hh"

namespace polytope {

class TriangleTessellator : public SerialTessellator<2> {
public:

  using RealType = double;
  using QuantizedTessellation = QuantTessellation<2>;
  using PrimitiveCells = VoronoiPrimitiveCells<2>;

  // Constructor, destructor.
  TriangleTessellator() = default;
  virtual ~TriangleTessellator() = default;

  // Compute the nodes around a collection of generators.
  // Required method for all Tessellators.
  PrimitiveCells tessellateQuantizedImpl(const QuantizedTessellation& result) const override;

  // The name of the tessellator
  virtual std::string name() const override { return "TriangleTessellator"; }

};

} //end polytope namespace

#endif
