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

  // Shifted coordinates in physical space
  RealPoint m_lx_o, m_xlo_o, m_dx_o;
  // Original coordinate locations
  RealPoint m_xlo, m_xhi;
  // Percent to pad the bounding box for good measure
  RealType m_pad = 0.5;
  // Maximum possible coordinate in a single direction
  // Not necessarily the max for this instance
  constexpr static IntType m_coordMax = HashKey<Dimension>::coordMax();
  // Current maximum coordinate in a single direction
  IntPoint maxCoord = IntPoint(m_coordMax);
  IntPoint minCoord = IntPoint();
  IntPoint maxBound = maxCoord - 1000;
  IntPoint minBound = minCoord + 1000;
  RealPoint rmaxBound = maxBound.template type_cast<RealType>();
  RealPoint rminBound = minBound.template type_cast<RealType>();
  bool m_init = false;

  Quantizer() = default;

  // Have the Quantizer determine the bounding box, for testing and serial mode
  // Combine point vectors if multiple exist
  Quantizer(const std::vector<RealType>& points,
            const RealType& degeneracy = -1.,
            const RealType& pad = -1.) {
    init(points, degeneracy, pad);
  }

  // Provide the min and max bounding points for the whole domain
  Quantizer(const RealPoint& xlo,
            const RealPoint& xhi,
            const RealType& degeneracy = -1.,
            const RealType& pad = -1.) {
    init(xlo, xhi, degeneracy, pad);
  }

  IntPoint quantize(const RealPoint& x) const {
    POLY_ASSERT2(m_init, "Must initialize quantizer before using it");
    return x.template convertXi<IntType>(m_xlo_o, m_dx_o);
  }

  RealPoint dequantize(const IntPoint& X) const {
    POLY_CHECK2(m_init, "Must initialize quantizer before using it");
    return X.convertx(m_xlo_o, m_dx_o);
  }

  CoordHash hash(const IntPoint& X) const {
    POLY_CHECK2(m_init, "Must initialize quantizer before using it");
    return Hasher::hash(X);
  }

  CoordHash hash_quantize(const RealPoint& x) const {
    POLY_CHECK2(m_init, "Must initialize quantizer before using it");
    return Hasher::hash(quantize(x));
  }

  IntPoint unhash(const CoordHash& h) const {
    POLY_CHECK2(m_init, "Must initialize quantizer before using it");
    return Hasher::unhash(h);
  }

  RealPoint unhash_dequantize(const CoordHash& h) const {
    POLY_CHECK2(m_init, "Must initialize quantizer before using it");
    return dequantize(unhash(h));
  }

  RealPoint degeneracy() const {
    return m_dx_o*m_lx_o;
  }

  void init(const RealPoint& xlo,
            const RealPoint& xhi,
            const RealType& degeneracy = -1.,
            const RealType& pad = -1.) {
    if (pad >= 0.) {
      m_pad = pad;
    }
    m_xlo = xlo;
    m_xhi = xhi;
    m_lx_o = (xhi - xlo)*(1.0 + m_pad);
    m_xlo_o = xlo - 0.5*(xhi - xlo)*m_pad;
    if (degeneracy > 0.) {
      m_dx_o = degeneracy*m_lx_o;
      maxCoord = (m_lx_o/m_dx_o).template type_cast<IntType>();
      maxBound = maxCoord - 1;
      rmaxBound = maxBound.template type_cast<RealType>();
    } else {
      m_dx_o = m_lx_o/static_cast<RealType>(m_coordMax);
    }
    m_init = true;
  }

  void init(const std::vector<RealType>& points,
            const RealType& degeneracy = -1.,
            const RealType& pad = -1.) {
    RealPoint minPoint(0.99*std::numeric_limits<RealType>::max());
    RealPoint maxPoint = -minPoint;
    std::vector<RealPoint> rpoints = extractCoords<Dimension, RealType>(points);
    findBoundingElements<Dimension, RealType>(rpoints, minPoint, maxPoint);
    init(minPoint, maxPoint, degeneracy, pad);
  }

  bool inBounds(const RealPoint& point) const {
    POLY_CHECK2(m_init, "Must initialize quantizer before using it");
    if (point.allLessEqual(m_xhi) && point.allGreaterEqual(m_xlo)) {
      return true;
    }
    return false;
  }

  bool inQBounds(const IntPoint& point) const {
    if (point.allLess(maxBound) && point.allGreater(minBound)) {
      return true;
    }
    return false;
  }

  bool inQBounds(const RealPoint& point) const {
    if (point.allLess(rmaxBound) && point.allGreater(rminBound)) {
      return true;
    }
    return false;
  }

  // Determine exactly which sides of a box a point is external to
  // 1:  Outside upper side
  // 0:  Not outside
  // -1: Outside lower side
  Point<Dimension, int> externalSides(const RealPoint& point) const {
    Point<Dimension, int> out;
    out.zero();
    for (int dir = 0; dir < Dimension; ++dir) {
      if (point[dir] > rmaxBound[dir] - 10.) {
        out[dir] = 1;
      } else if (point[dir] < rminBound[dir] + 10.) {
        out[dir] = -1;
      }
    }
    return out;
  }
};
} // namespace polytope
#endif
