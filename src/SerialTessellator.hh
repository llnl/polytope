//----------------------------------------------------------------------------//
// SerialTessellator
//
// Common serial tessellator implementation.  Backends provide raw Voronoi
// primitives; VoronoiAssembler performs the shared clipping and topology
// construction needed to produce a QuantTessellation.
//----------------------------------------------------------------------------//
#ifndef __Polytope_SerialTessellator__
#define __Polytope_SerialTessellator__

#include "Tessellator.hh"
#include "VoronoiAssembler.hh"

namespace polytope {

template<int Dimension>
class SerialTessellator : public Tessellator<Dimension, double> {
public:
  using Base = Tessellator<Dimension, double>;
  using QuantizedTessellation = QuantTessellation<Dimension>;
  using PrimitiveCells = VoronoiPrimitiveCells<Dimension>;

  virtual ~SerialTessellator() = default;

  void tessellateQuantized(QuantizedTessellation& result) final override {
    if (result.points.empty()) {
      return;
    } else if (result.points.size() == 1) {
      this->singleNodeTessellate(result);
    } else {
      VoronoiAssembler<Dimension> assembler(result);
      assembler.assemble(this->tessellateQuantizedImpl(result));
    }
  }

protected:
  //! Generate one collection of raw Voronoi primitives per generator.
  virtual PrimitiveCells
  tessellateQuantizedImpl(const QuantizedTessellation& input) const = 0;
};

} // namespace polytope

#endif
