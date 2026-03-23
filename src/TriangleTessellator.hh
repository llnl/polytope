//------------------------------------------------------------------------
// TriangleTessellator
// 
// An implemenation of the Tessellator interface that uses the Triangle
// library by Jonathan Shewchuk.
//------------------------------------------------------------------------
#ifndef __Polytope_TriangleTessellator__
#define __Polytope_TriangleTessellator__

#ifdef POLYTOPE_ENABLE_TRIANGLE

#include <vector>
#include <cmath>

#include "Tessellator.hh"
#include "QuantizedTessellation2d.hh"
#include "Point.hh"
#include "polytope_tessellator_utilities.hh"

namespace polytope {

template<typename RealType>
class TriangleTessellator: public Tessellator<2, RealType> {
public:

  // Typedefs for edges, coordinates, and points
  typedef int                 CoordHash;
  typedef std::pair<int, int> EdgeHash;
  typedef typename Tessellator<2, RealType>::QuantizedTessellation QuantizedTessellation;
  
  // The typedefs that follow from this choice
  typedef Point2<RealType>  RealPoint;
  typedef Point2<CoordHash> IntPoint;

  // Constructor, destructor.
  TriangleTessellator();
  ~TriangleTessellator();

  // Compute the nodes around a collection of generators.
  // Required method for all Tessellators.
  virtual void tessellateQuantized(QuantizedTessellation& result) const;

  // Return the tessellator name
  std::string name() const { return "TriangleTessellator"; }

  //! Returns the accuracy to which this tessellator can distinguish coordinates.
  //! Should be returned appropriately for normalized coordinates, i.e., if all
  //! coordinates are in the range xi \in [0,1], what is the minimum allowed 
  //! delta in x.
  virtual RealType degeneracy() const { return mDegeneracy; }
  void degeneracy(RealType degeneracy) const { mDegeneracy = degeneracy; }

private:
  //-------------------- Private interface ---------------------- //
  static RealType mDegeneracy; 
};


//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
template<typename RealType> 
RealType  
TriangleTessellator<RealType>::mDegeneracy = 8.0/std::numeric_limits<typename TriangleTessellator<RealType>::CoordHash>::max();

} //end polytope namespace

#endif
#endif
