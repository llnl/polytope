//----------------------------------------------------------------------------//
// SerialTessellator
//
// Common serial tessellator implementation. Backends provide a Voronoi
// assembler; which is then used to perform the shared clipping and topology
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

  virtual ~SerialTessellator() = default;

  void tessellateQuantized(QuantizedTessellation& result) final override {
    if (result.points.empty()) {
      return;
    } else if (result.points.size() == 1) {
      this->singleNodeTessellate(result);
    } else {
      VoronoiAssembler<Dimension> assembler = this->tessellateQuantizedImpl(result);
      assembler.assemble();
    }
  }

  //! Generate one collection of raw Voronoi primitives per generator.
  virtual VoronoiAssembler<Dimension>
  tessellateQuantizedImpl(QuantizedTessellation& input) const = 0;
};

} // namespace polytope

#endif
