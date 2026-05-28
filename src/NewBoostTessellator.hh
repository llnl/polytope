//------------------------------------------------------------------------
// BoostTessellator
// 
// Polytope wrapper for the native 2D Voronoi tessellator in Boost.Polygon
// v1.52 or greater
//------------------------------------------------------------------------
#ifndef __Polytope_BoostTessellator__
#define __Polytope_BoostTessellator__

#include <vector>
#include <cmath>
#include <limits>

#include "boost/polygon/voronoi.hpp"

#include "Tessellator.hh"
#include "QuantTessellation.hh"
#include "Point.hh"
//#include "polytope_tessellator_utilities.hh"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/algorithms/unique.hpp>
#include <boost/bind/bind.hpp>
//BOOST_GEOMETRY_REGISTER_POINT_2D(polytope::Point2<double>, double, boost::geometry::cs::cartesian, x, y);
// BOOST_GEOMETRY_REGISTER_POINT_2D(polytope::Point2<int32_t>, int32_t, boost::geometry::cs::cartesian, x, y);
// BOOST_GEOMETRY_REGISTER_POINT_2D(polytope::Point2<int64_t>, int64_t, boost::geometry::cs::cartesian, x, y);
namespace polytope {

template<typename RealType>
class BoostTessellator: public Tessellator<2, RealType> {
public:

  // The Boost.Polygon Voronoi diagram
  typedef boost::polygon::voronoi_diagram<RealType> VD;
  using QuantizedTessellation = QuantTessellation<2>;

  // Constructor, destructor.
  BoostTessellator() = default;
  virtual ~BoostTessellator() = default;

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
  RealType mDegeneracy = 0.;
};

//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
// template<typename RealType> 
// RealType  
// BoostTessellator<RealType>::mDegeneracy = 8.0/std::numeric_limits<typename BoostTessellator<RealType>::CoordHash>::max();

} //end polytope namespace

#endif
