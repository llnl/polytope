//------------------------------------------------------------------------
// TetgenTessellator
// 
// An implemenation of the Tessellator interface that uses the Tetgen
// library by Hang Si.
// By default tetgen is built assuming double coordinates, so we only
// provide that implementation as our RealType.
//------------------------------------------------------------------------
#ifndef __Polytope_TetgenTessellator__
#define __Polytope_TetgenTessellator__

#ifdef POLYTOPE_ENABLE_TETGEN

#include <vector>
#include <cmath>

#include "Tessellator.hh"
#include "Point.hh"
#include "QuantTessellation.hh"
#include "ReducedPLC.hh"

namespace polytope {

class TetgenTessellator: public Tessellator<3, double> {
public:
  //-------------------- Public interface --------------------
  typedef double RealType;
  typedef Tessellator<3, RealType>::QuantizedTessellation QuantizedTessellation;

  // Constructor, destructor.
  TetgenTessellator();
  ~TetgenTessellator();

  // Tessellate the given generators (unbounded).
  void tessellate(const std::vector<RealType>& points,
                  Tessellation<3, RealType>& mesh) const;

  // Tessellate with a bounding box representing the boundaries.
  void tessellate(const std::vector<RealType>& points,
                  RealType* low,
                  RealType* high,
                  Tessellation<3, RealType>& mesh) const;

  // Tessellate obeying the given boundaries.
  void tessellate(const std::vector<RealType>& points,
                  const std::vector<RealType>& PLCpoints,
                  const PLC<3, RealType>& geometry,
                  Tessellation<3, RealType>& mesh) const;

  // Tessellate obeying the given boundaries expressed as a ReducedPLC.
  void tessellate(const std::vector<RealType>& points,
                  const ReducedPLC<3, RealType>& geometry,
                  Tessellation<3, RealType>& mesh) const;

  // This Tessellator handles PLCs!
  bool handlesPLCs() const { return true; }

  // This backend still uses the older direct tessellation path in 3D.
  // We define the method so the class is concrete, but the quantized API
  // migration remains separate work.
  void tessellateQuantized(QuantizedTessellation& qmesh) const;

  // Return the name of this tessellator
  std::string name() const { return "TetgenTessellator"; }

  //! Returns the accuracy to which this tessellator can distinguish coordinates.
  //! Should be returned appropriately for normalized coordinates, i.e., if all
  //! coordinates are in the range xi \in [0,1], what is the minimum allowed 
  //! delta in x.
  virtual RealType degeneracy() const { return mDegeneracy; }

private:
  //-------------------- Private interface --------------------
  typedef int64_t CoordHash;
  static CoordHash coordMax;
  static RealType mDegeneracy;

  // Internal method to compute the quantized tessellation.
  void computeUnboundedQuantizedTessellation(const std::vector<RealType>& points,
                                             const std::vector<RealType>& nonGeneratingPoints,
                                             internal::QuantTessellation<3, RealType>& mesh) const;
};

}

#endif
#endif
