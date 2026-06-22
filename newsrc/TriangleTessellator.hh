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

#include "Tessellator.hh"
#include "QuantTessellation.hh"
#include "Point.hh"

namespace polytope {

class TriangleTessellator : public Tessellator<2, double> {
public:

  using RealType = double;
  using QuantizedTessellation = QuantTessellation<2>;

  // Constructor, destructor.
  TriangleTessellator() = default;
  TriangleTessellator(const Quantizer<2>& Q) :
    Tessellator(Q) {}
  virtual ~TriangleTessellator() = default;

  // Compute the nodes around a collection of generators.
  // Required method for all Tessellators.
  virtual void tessellateQuantized(const QuantPLC<2>& qplc,
                                   QuantizedTessellation& result) const;

  // The name of the tessellator
  std::string name() const { return "TriangleTessellator"; }

};

} //end polytope namespace

#endif
