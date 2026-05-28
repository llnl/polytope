//-----------------------------------------------------------------------------//
// QuantTessellation
//
//-----------------------------------------------------------------------------//

#ifndef POLYTOPE_QUANTTESSELLATION_HH
#define POLYTOPE_QUANTTESSELLATION_HH

#include <vector>
#include "HashKey.hh"
#include "Point.hh"
#include "Quantizer.hh"
#include "Tessellation.hh"

namespace polytope {

template<int Dimension>
class QuantTessellation {
public:
  using RealType = double;
  using CoordHash = typename HashKey<Dimension>::CoordHash;
  using IntType = typename HashKey<Dimension>::IntType;
  using IntPoint = Point<Dimension, IntType>;
  using RealPoint = Point<Dimension, RealType>;
  using Quant = Quantizer<Dimension>;
  using Hasher = HashKey<Dimension>;
  using TessellationType = Tessellation<Dimension, RealType>;

  QuantTessellation(const Quant& Q,
                    const std::vector<RealType>& genpoints);

  // Returns quantized points cast as reals to give to the tessellator
  std::vector<RealPoint> getRealPoints() const;
  // Returns quantized points cast as reals to give to the tessellator
  std::vector<IntPoint> getIntPoints() const;

  // Dequantize and fill the tessellation type after tessellation
  void fillTessellation(TessellationType& mesh) const;

  //------------------------------------------------------------------------------
  // Member data
  //------------------------------------------------------------------------------
  const Quant m_Q;
  // Generator points
  std::vector<CoordHash> m_hashes;
  std::vector<IntPoint> m_points;
  // Local lower and upper bounding box coordinates
  IntPoint m_loBounds;
  IntPoint m_hiBounds;
  std::vector<IntPoint> m_nodes; // Nodes that make up the Voronoi
  std::vector<std::vector<int>> m_faces; // Faces made up of indices into m_nodes
  std::vector<std::vector<int>> m_cells; // Cells made up of indices into m_faces
private:
  template<int D = Dimension>
  std::enable_if_t<D == 2, void>
  fillTessellation2D(TessellationType& mesh) const;
  template<int D = Dimension>
  std::enable_if_t<D == 3, void>
  fillTessellation3D(TessellationType& mesh) const;
};

} // namespace polytope
#endif
