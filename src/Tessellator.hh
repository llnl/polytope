#ifndef __Polytope_Tessellator__
#define __Polytope_Tessellator__

#include "QuantTessellation.hh"
#include "polytope_internal.hh"
#include "Tessellation.hh"

#include <optional>

namespace polytope {

//! \class Tessellator - An abstract base class for objects that generate
//! Voronoi and Voronoi-like tessellations for sets of points and/or
//! geometries.
template<int Dimension, typename RealType>
class Tessellator {
public:
  using TessType = Tessellation<Dimension, RealType>;
  using RealPoint = Point<Dimension, RealType>;
  using RealPointVec = std::vector<RealPoint>;

  //! Default constructor.
  Tessellator() = default;

  //! Destructor.
  virtual ~Tessellator() {}

  //! Generate a Voronoi tessellation for the given set of generator points.
  //! \param points A (numPoints) array containing point coordinates.
  //! \return A self-contained, real-coordinate tessellation.
  virtual TessType tessellate(const std::vector<Point<Dimension, RealType>>& points);

  //! Wrapper for above virtual function
  TessType tessellate(const std::vector<RealType>& points) {
    return this->tessellate(extractCoords<Dimension, RealType>(points));
  }

  //! Generate a Voronoi-like tessellation for the given set of generator
  //! points and a description of the geometry in which they exist.
  //! \param points A (numPoints) array containing point coordinates.
  //! \param PLCpoints A (n) array containing point coordinates for the PLC.
  //! \param geometry A description of the geometry in Piecewise Linear Complex form.
  //! \return A self-contained, real-coordinate tessellation.
  virtual TessType tessellate(const std::vector<Point<Dimension, RealType>>& points,
                              const std::vector<Point<Dimension, RealType>>& PLCpoints,
                              const PLC<Dimension>& geometry);

  //! Wrapper for above virtual function
  TessType tessellate(const std::vector<RealType>& points,
                      const std::vector<RealType>& PLCpoints,
                      const PLC<Dimension>& geometry) {
    return this->tessellate(extractCoords<Dimension, RealType>(points),
                            extractCoords<Dimension, RealType>(PLCpoints),
                            geometry);
  }

  //! Wrapper for above virtual function
  TessType tessellate(const std::vector<Point<Dimension, RealType>>& points,
                      const std::vector<RealType>& PLCpoints,
                      const PLC<Dimension>& geometry) {
    return this->tessellate(points, extractCoords<Dimension, RealType>(PLCpoints),
                            geometry);
  }

  //! Compute the quantized tessellation.  Serial implementations inherit the
  //! primitive-producing implementation from SerialTessellator; distributed
  //! implementations provide this operation directly.
  virtual void
  tessellateQuantized(QuantTessellation<Dimension>& result) = 0;

  //! Required for all tessellators:
  //! A unique name string per tessellation instance.
  virtual std::string name() const = 0;

  //! Required for all tessellators:
  //! Returns the accuracy to which this tessellator can distinguish coordinates.
  //! Should be returned appropriately for normalized coordinates, i.e., if all
  //! coordinates are in the range xi \in [0,1], what is the minimum allowed
  //! delta in x.
  Point<Dimension, double> degeneracy() const {
    return m_Q.degeneracy();
  }

  //! Set the padding for quantized space. Reinitializes quantizer if necessary
  void setPadding(const RealType& pad) {
    m_Q.setPadding(pad);
  }

  RealType getPadding() const {
    return m_Q.getPadding();
  }

  Quantizer<Dimension>& quantizer() {
    return m_Q;
  }

  const Quantizer<Dimension>& quantizer() const {
    return m_Q;
  }

  void useMortonEncoding() { m_Q.useMortonEncoding(); }
  void usePackedEncoding() { m_Q.usePackedEncoding(); }

  void singleNodeTessellate(QuantTessellation<Dimension>& quantmesh);

  bool clippingEnabled() const {
#ifdef POLYTOPE_ENABLE_HIBIT2D
    return false;
#else
    return true;
#endif
  }

  void setQBounds(const RealPointVec& qbounds) {
    m_QBounds.emplace(qbounds);
  }

protected:
  Quantizer<Dimension> m_Q;

  // Initialize the Q. Check if QBounds has been provided.
  void initQ(RealPointVec points,
             const RealPointVec& PLCpoints) {
    points.insert(points.end(), PLCpoints.begin(), PLCpoints.end());
    initQ(points);
  }

  void initQ(RealPointVec points) {
    if (m_QBounds.has_value()) {
      points.insert(points.end(), m_QBounds->begin(), m_QBounds->end());
    }
    m_Q.init(points);
  }

private:
  // Optional set of bounding points to give to quantizer
  // These are not used as generators
  std::optional<RealPointVec> m_QBounds;
};

}

#include "TessellatorInline.hh"

#endif
