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
//    is considered visible if it's Voronoi cell intersects the convex hull at all.
// 7. The entire set of visible generators are communicated to all processors.
// 8. Each proc creates a quantized Voronoi using all visible generators (IE the visible Voronoi).
//    This is same on all processors. Each generator retains which processor it originated from.
// 9. Each proc determines which processors are neighbors. Specifically, procs A and B
//    are neighbors if either is true:
//    a. If the convex hull of proc A intersects with the convex hull of proc B.
//    b. In the visible Voronoi, if a cell created from a generator from proc A shares
//       a face with a cell created from a generator from proc B.
// 10. Each proc exchanges the ENTIRE set of generators with any neighbors (not ideal).
//     So proc A would get all generators from proc B.
// 11. Each proc generates a Voronoi diagram using it's generators and neighbor generators,
//     but not the visible generators.
// 12. Clip the Voronoi generated from both local and neighbor generators.
// 13. Filter out any points and cells that are not local to that rank.
//----------------------------------------------------------------------------//
#ifndef __Polytope_DistributedTessellator__
#define __Polytope_DistributedTessellator__

#include <string>

#include "polytope.hh"
#include "Partitioner.hh"
#include "Tessellator.hh"

namespace polytope {

template<int Dimension>
class DistributedTessellator: public Tessellator<Dimension, double> {

  //--------------------------- Public Interface ---------------------------//
public:
  using RealType = double;
  using Base = Tessellator<Dimension, RealType>;
  using QuantizedTessellation = QuantTessellation<Dimension>;
  using TessellationType = Tessellation<Dimension, RealType>;

  DistributedTessellator(Base& serialTessellator);
  virtual ~DistributedTessellator() {
    m_QPLC.clear();
    m_clipping = false;
  }

  //! Select whether generator exchanges send quantized points instead of keys.
  void setExchangePoints(const bool value) { m_exchangePoints = value; }
  bool exchangePoints() const { return m_exchangePoints; }

  virtual void tessellate(const std::vector<RealType>& points,
                          TessellationType& mesh) override;

  virtual void tessellate(const std::vector<RealType>& points,
                          const std::vector<RealType>& PLCpoints,
                          const PLC<Dimension>& geometry,
                          TessellationType& mesh) override;

  //! Partition a replicated generator set, then tessellate the local result.
  //! Every rank must pass the same points in the same order.
  void partitionAndTessellate(const std::vector<RealType>& points,
                              const Partitioner<Dimension>& partitioner,
                              TessellationType& mesh);

  //! Partition a replicated generator set inside a bounding PLC, then
  //! tessellate the local result. Every rank must pass the same points in the
  //! same order.
  void partitionAndTessellate(const std::vector<RealType>& points,
                              const std::vector<RealType>& PLCpoints,
                              const PLC<Dimension>& geometry,
                              const Partitioner<Dimension>& partitioner,
                              TessellationType& mesh);

  //! Simply becomes a wrapper for the Impl
  virtual void tessellateQuantized(QuantizedTessellation& qmesh) override {
    this->tessellateQuantizedImpl(qmesh);
  }

  virtual void tessellateQuantizedImpl(QuantizedTessellation& qmesh) override;

  virtual std::string name() const override;

  QuantizedTessellation generateVisibleMesh(QuantizedTessellation& qmesh);

  void checkEncoding() {
    auto qkey = Quantizer<Dimension>::instance().keyEncoding();
    POLY_CHECK2(qkey == m_keyEncode,
                "Key encoding method has changed during tessellation");
  }

private:
  //! Broadcast the root rank's requested exchange representation.
  void synchronizeExchangePoints();

  Base& m_serialTessellator;
  // Whether to exchange points in parallel or exchange keys
  bool m_exchangePoints = false;
  KeyEncoding m_keyEncode;
  bool m_clipping = false;
  QuantPLC<Dimension> m_QPLC;

  // Forbidden methods.
  DistributedTessellator();
  DistributedTessellator(const DistributedTessellator&);
  DistributedTessellator& operator=(const DistributedTessellator&);
};

}

#endif
