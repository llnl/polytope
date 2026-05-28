//-----------------------------------------------------------------------------//
// Quantizer
//
// Class for quantizing and dequantizing points or PLCs.
//-----------------------------------------------------------------------------//
#ifndef __Polytope_Quantizer__
#define __Polytope_Quantizer__

#include "Point.hh"
#include "HashKey.hh"

namespace polytope {

template<int Dimension>
class Quantizer {
public:
  using RealType = double;
  using Hasher = HashKey<Dimension>;
  using CoordHash = typename Hasher::CoordHash;
  using IntType = typename Hasher::IntType;
  using IntPoint = Point<Dimension, IntType>;
  using RealPoint = Point<Dimension, RealType>;

  // Lower case represents coordinates in physical space
  // Upper case represents coordinates in transformed space
  // Subscript i means inner box and o means outerbox
  RealPoint m_lx_o, m_xlo_o, m_dx_o;
  // Percent to pad the bounding box for good measure
  RealType m_pad = 0.02;
  // TODO: Implement smaller box logic
  // IntPoint m_Xlo_i, m_Xhi_i;
  // IntType m_dX_i = 0;
  constexpr static IntType coordMax = HashKey<Dimension>::coordMax();

  // Provide the min and max bounding points for the whole domain
  Quantizer(RealPoint& xlo, RealPoint& xhi, RealType pad = 0.02) :
    m_lx_o((xhi - xlo)*(1.0 + pad)),
    m_xlo_o(xlo - 0.5*(xhi - xlo)*pad),
    m_dx_o(m_lx_o/static_cast<RealType>(coordMax)),
    m_pad(pad) { }

  // Provide the bounding points for the inner and outer boxes
  // Quantizer(RealPoint& xlo, RealPoint& xhi,
  //           RealPoint& xlo_i, RealPoint& xhi_i) :
  //   Quantizer(xlo, xhi) { }

  IntPoint quantize(const RealPoint& x) const {
    return x.template convertXi<IntType>(m_xlo_o, m_dx_o);
  }

  RealPoint dequantize(const IntPoint& X) const {
    return X.convertx(m_xlo_o, m_dx_o);
  }

  CoordHash hash(const IntPoint& X) const {
    return Hasher::hash(X);
  }

  CoordHash hash_quantize(const RealPoint& x) const {
    return Hasher::hash(quantize(x));
  }

  IntPoint unhash(const CoordHash& h) const {
    return Hasher::unhash(h);
  }

  RealPoint unhash_dequantize(const CoordHash& h) const {
    return dequantize(unhash(h));
  }

  RealPoint degeneracy() const {
    return m_dx_o;
  }

};
} // namespace polytope
#endif
