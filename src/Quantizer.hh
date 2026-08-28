//-----------------------------------------------------------------------------//
// Quantizer
//
// Singleton class for quantizing and dequantizing points or PLCs.
//-----------------------------------------------------------------------------//
#ifndef __Polytope_Quantizer__
#define __Polytope_Quantizer__

#include "Point.hh"
#include "KeyCodec.hh"
#include "Communicator.hh"
#include <mutex>

namespace polytope {

template<int Dimension>
class Quantizer {
public:
  using RealType = double;
  using Traits = QuantizedKeyTraits<Dimension>;
  using Codec = KeyCodec<Dimension>;
  using Key = typename Codec::Key;
  using RealPoint = Point<Dimension, RealType>;

  // Delete copy constructor and assignment operator
  Quantizer(const Quantizer&) = delete;
  Quantizer& operator=(const Quantizer&) = delete;

  // Get the singleton instance
  static Quantizer& instance() {
    static Quantizer instance;
    return instance;
  }

  KeyEncoding keyEncoding() const {
    return m_codec.encoding();
  }

  const std::string keyName() const noexcept {
    return m_codec.keyName();
  }

  // Shifted coordinates in physical space
  RealPoint m_lx_o, m_xlo_o, m_dx_o;
  // Original coordinate locations
  RealPoint m_xlo, m_xhi;
  // Percent to pad the bounding box for good measure
  RealType m_pad = 0.04;
  // Maximum possible coordinate in a single direction
  // Not necessarily the max for this instance
  constexpr static QuantizedCoordinate<Dimension> m_maxCoordinate = Traits::maxCoordinate();
  // Current maximum coordinate in a single direction
  QuantizedPoint<Dimension> maxCoord = QuantizedPoint<Dimension>(m_maxCoordinate);
  QuantizedPoint<Dimension> minCoord = QuantizedPoint<Dimension>();
  QuantizedPoint<Dimension> maxBound = maxCoord - 1000;
  QuantizedPoint<Dimension> minBound = minCoord + 1000;
  RealPoint rmaxBound = maxBound.template type_cast<RealType>();
  RealPoint rminBound = minBound.template type_cast<RealType>();
  bool m_init = false;

  // Extend or reduce the padding on the boundary box by a certain percentage
  void extend(const RealType extendPad) {
    std::lock_guard<std::mutex> lock(m_mutex);
    POLY_CHECK2(m_init, "Must initialize quantizer before extending it");
    _init_impl(m_xlo, m_xhi, -1, extendPad);
  }

  QuantizedPoint<Dimension> quantize(const RealPoint& x) const {
    POLY_ASSERT2(m_init, "Must initialize quantizer before using it");
    return x.template convertXi<QuantizedCoordinate<Dimension>>(m_xlo_o, m_dx_o);
  }

  RealPoint dequantize(const QuantizedPoint<Dimension>& X) const {
    POLY_CHECK2(m_init, "Must initialize quantizer before using it");
    return X.convertx(m_xlo_o, m_dx_o);
  }

  Key encode(const QuantizedPoint<Dimension>& X) const {
    POLY_CHECK2(m_init, "Must initialize quantizer before using it");
    return m_codec.encode(X);
  }

  Key quantizeAndEncode(const RealPoint& x) const {
    POLY_CHECK2(m_init, "Must initialize quantizer before using it");
    return m_codec.encode(quantize(x));
  }

  QuantizedPoint<Dimension> decode(const Key& h) const {
    POLY_CHECK2(m_init, "Must initialize quantizer before using it");
    return m_codec.decode(h);
  }

  RealPoint decodeAndDequantize(const Key& h) const {
    POLY_CHECK2(m_init, "Must initialize quantizer before using it");
    return dequantize(decode(h));
  }

  RealPoint degeneracy() const {
    return m_dx_o*m_lx_o;
  }

  void useMortonEncoding() {
    m_codec.setEncoding(KeyEncoding::Morton);
  }

  void usePackedEncoding() {
    m_codec.setEncoding(KeyEncoding::Packed);
  }

  void init(const RealPoint& xlo,
            const RealPoint& xhi,
            const RealType& degeneracy = -1.,
            const RealType& pad = -1.) {
    std::lock_guard<std::mutex> lock(m_mutex);
    _init_impl(xlo, xhi, degeneracy, pad);
  }

  void init(const std::vector<RealType>& points,
            const RealType& degeneracy = -1.,
            const RealType& pad = -1.) {
    std::lock_guard<std::mutex> lock(m_mutex);
    RealPoint minPoint(0.99*std::numeric_limits<RealType>::max());
    RealPoint maxPoint = -minPoint;
    std::vector<RealPoint> rpoints = extractCoords<Dimension, RealType>(points);
    findBoundingElements<Dimension, RealType>(rpoints, minPoint, maxPoint);
    _init_impl(minPoint, maxPoint, degeneracy, pad);
  }

  bool inBounds(const RealPoint& point) const {
    POLY_CHECK2(m_init, "Must initialize quantizer before using it");
    if (point.allLessEqual(m_xhi) && point.allGreaterEqual(m_xlo)) {
      return true;
    }
    return false;
  }

  bool inQBounds(const QuantizedPoint<Dimension>& point) const {
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

private:
  // Private constructor
  Quantizer() = default;

  // Private initialization implementation
  void _init_impl(const RealPoint& xlo,
                  const RealPoint& xhi,
                  const RealType degeneracy,
                  const RealType pad) {
    if (pad >= 0.) {
      m_pad = pad;
    }
    m_xlo = xlo;
    m_xhi = xhi;
#ifdef POLYTOPE_ENABLE_MPI
    auto comm = Communicator::communicator();
    auto localMin = xlo.toArray();
    auto localMax = xhi.toArray();
    std::array<RealType, Dimension> reducedMin, reducedMax;
    MPI_Allreduce(localMin.data(), reducedMin.data(), Dimension, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(localMax.data(), reducedMax.data(), Dimension, MPI_DOUBLE, MPI_MAX, comm);
    m_xlo = RealPoint(reducedMin.data());
    m_xhi = RealPoint(reducedMax.data());
#endif
    m_lx_o = (xhi - xlo)*(1.0 + m_pad);
    m_xlo_o = xlo - 0.5*(xhi - xlo)*m_pad;
    if (degeneracy > 0.) {
      m_dx_o = degeneracy*m_lx_o;
      maxCoord = (m_lx_o/m_dx_o).template type_cast<QuantizedCoordinate<Dimension>>();
      maxBound = maxCoord - 1;
      rmaxBound = maxBound.template type_cast<RealType>();
    } else {
      m_dx_o = m_lx_o/static_cast<RealType>(m_maxCoordinate);
    }
    m_init = true;
  }

  // Mutex for thread-safe initialization
  Codec m_codec;
  mutable std::mutex m_mutex;
};

} // namespace polytope
#endif
