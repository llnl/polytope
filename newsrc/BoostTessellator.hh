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

#include "Tessellator.hh"
#include "QuantTessellation.hh"
#include "Point.hh"

namespace polytope {

class BoostTessellator : public Tessellator<2, double> {
public:

  using RealType = double;
  using QuantizedTessellation = QuantTessellation<2>;

  // Constructor, destructor.
  BoostTessellator() = default;
  virtual ~BoostTessellator() = default;

  // Compute the nodes around a collection of generators.
  // Required method for all Tessellators.
  virtual void tessellateQuantized(QuantizedTessellation& result) const;

  // The name of the tessellator
  std::string name() const { return "BoostTessellator"; }

};

} //end polytope namespace

#endif
