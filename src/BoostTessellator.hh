//------------------------------------------------------------------------
// BoostTessellator
// 
// Polytope wrapper for the native 2D Voronoi tessellator in Boost.Polygon
// v1.52 or greater
//------------------------------------------------------------------------
#ifndef __Polytope_BoostTessellator__
#define __Polytope_BoostTessellator__

#ifdef POLYTOPE_ENABLE_BOOST

#include <vector>
#include <cmath>
#include <limits>

#include "boost/polygon/voronoi.hpp"

#include "Tessellator.hh"
#include "QuantizedTessellation2d.hh"
#include "Point.hh"
#include "polytope_tessellator_utilities.hh"

namespace polytope {

template<typename RealType>
class BoostTessellator: public Tessellator<2, RealType> {
public:

  // The Boost.Polygon Voronoi diagram
  typedef boost::polygon::voronoi_diagram<RealType> VD;

  // Some useful typedefs
  typedef int                 CoordHash;
  typedef std::pair<int, int> EdgeHash;
  typedef typename Tessellator<2, RealType>::QuantizedTessellation QuantizedTessellation;

  // The typedefs that follow from this choice
  typedef Point2<RealType>  RealPoint;
  typedef Point2<CoordHash> IntPoint;

  // Constructor, destructor.
  BoostTessellator();
  ~BoostTessellator();

  // Compute the nodes around a collection of generators.
  // Required method for all Tessellators.
  virtual void tessellateQuantized(QuantizedTessellation& result) const;

  // The name of the tessellator
  std::string name() const { return "BoostTessellator"; }

  //! Returns the accuracy to which this tessellator can distinguish coordinates.
  //! Should be returned appropriately for normalized coordinates, i.e., if all
  //! coordinates are in the range xi \in [0,1], what is the minimum allowed 
  //! delta in x.
  virtual RealType degeneracy() const { return mDegeneracy; }
  void degeneracy(const RealType val) const { mDegeneracy = val; }

private:
  //-------------------- Private interface ---------------------- //
  static RealType mDegeneracy; 
};

//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
template<typename RealType> 
RealType  
BoostTessellator<RealType>::mDegeneracy = 8.0/std::numeric_limits<typename BoostTessellator<RealType>::CoordHash>::max();

} //end polytope namespace

#endif
#endif
