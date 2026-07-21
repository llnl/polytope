//----------------------------------------------------------------------------//
// DistributedTessellator
//
// Provides a parallel tessellation.
//
// Based on the parallel tessellation algorithm originally implemented in
// Spheral++.
// Steps for doing parallel tessellation:
// 1. Determine the boundary extents of either the generator points or the bounding PLC
// 2. Exchange these extents and initialize the Quantizer on each proc.
// 3. Create a QuantTessellation on each processor. All computations should be done
//    in quantized space as much as possible.
// 3*. Possibly redistribute generators so each local set is spatially contained.
// 4. Each processor creates a quantized voronoi of it's set of local generators (g_p).
// 5. The g_p is used to create a convex hull (as a QuantPLC) for each proc.
// 6. Each proc determines a set of local visible generators (v_p). A generator
//    is considered visible if it's voronoi cell intersects the convex at all.
// 7. The entire set of visible generators are communicated to all processors.
// 8. Each proc creates a quantized voronoi using all visible generators.
// 9. Each proc determines which set of generators are neighbors if either are true:
//    a. Does the convex hull overlap with another convex hull?
//    b. Are any cells between procs neighbors in the visible generator voronoi?
// 10. Each proc exchanges the ENTIRE set of generators with any neighbors (not ideal).
// 11. Each proc generates a voronoi diagram using it's generators and neighbor generators.
//----------------------------------------------------------------------------//
#ifndef __Polytope_DistributedTessellator__
#define __Polytope_DistributedTessellator__

#include <string>

#include "polytope.hh"
#include "Tessellator.hh"
#include "Communicator.hh"

namespace polytope {

template<int Dimension>
class DistributedTessellator: public Tessellator<Dimension, double> {

  //--------------------------- Public Interface ---------------------------//
public:
  using RealType = double;
  using Base = Tessellator<Dimension, RealType>;
  using QuantizedTessellation = QuantTessellation<Dimension>;
  using TessellationType = Tessellation<Dimension, RealType>;

#ifdef POLYTOPE_ENABLE_MPI
  DistributedTessellator(const Base& serialTessellator);
  virtual ~DistributedTessellator() = default;

  virtual void tessellate(const std::vector<RealType>& points,
                          TessellationType& mesh) const override;

  virtual void tessellate(const std::vector<RealType>& points,
                          const std::vector<RealType>& PLCpoints,
                          const PLC<Dimension>& geometry,
                          TessellationType& mesh) const override;

  //! Simply becomes a wrapper for the Impl
  virtual void tessellateQuantized(QuantizedTessellation& result) const override {
    this->tessellateQuantizedImpl(result);
  }
  virtual void tessellateQuantizedImpl(QuantizedTessellation& result) const override;

  virtual std::string name() const override;
#endif

private:
  const Base& m_serialTessellator;

  // Forbidden methods.
  DistributedTessellator();
  DistributedTessellator(const DistributedTessellator&);
  DistributedTessellator& operator=(const DistributedTessellator&);
};

}

#endif
