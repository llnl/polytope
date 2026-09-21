//------------------------------------------------------------------------
// BoostTessellator
// 
// Polytope wrapper for the native 2D Voronoi tessellator in Boost.Polygon
// v1.52 or greater
//------------------------------------------------------------------------
#ifndef __Polytope_BoostTessellator__
#define __Polytope_BoostTessellator__

#include <vector>
#include <cmath>
#include <limits>

#include "boost/polygon/voronoi.hpp"

#include "SerialTessellator.hh"
#include "QuantTessellation.hh"
#include "Point.hh"
#include "VoronoiAssembler.hh"

namespace polytope {

class BoostTessellator : public SerialTessellator<2> {
public:

  using RealType = double;
  using QuantizedTessellation = QuantTessellation<2>;
  using PrimitiveCells = VoronoiPrimitiveCells<2>;

  // Constructor, destructor.
  BoostTessellator() = default;
  virtual ~BoostTessellator() = default;

  // Compute the nodes around a collection of generators.
  // Required method for all Tessellators.
  PrimitiveCells tessellateQuantizedImpl(const QuantizedTessellation& result) const override;

  // The name of the tessellator
  virtual std::string name() const override { return "BoostTessellator"; }

};

} //end polytope namespace

#endif
