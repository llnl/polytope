#ifndef POLYTOPE_BOOST_UTILITIES_HH
#define POLYTOPE_BOOST_UTILITIES_HH
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include "PLC.hh"
//------------------------------------------------------------------------------
// A collection of utilities for doing things with Boost
//------------------------------------------------------------------------------

namespace polytope {

// We use the Boost.Geometry library to handle polygon intersections and such.
namespace bg = boost::geometry;

template <typename RealType, int Dimension>
using BGPoint = bg::model::point<RealType, Dimension, bg::cs::cartesian>;

template <typename RealType, int Dimension>
using BGPolygon = bg::model::polygon<BGPoint<RealType, Dimension>,false>;

//------------------------------------------------------------------------------
// Make a 2D Boost.Geometry point
//------------------------------------------------------------------------------
template <typename RealType>
BGPoint<RealType,2> makeBGPoint2D(std::vector<RealType>& position) {
   POLY_CHECK( position.size() == 2 );
   BGPoint<RealType,2> point(position[0], position[1]);
   return point;
}

//------------------------------------------------------------------------------
// Make a 3D Boost.Geometry point
//------------------------------------------------------------------------------
template <typename RealType>
BGPoint<RealType,3> makeBGPoint3D(std::vector<RealType>& position) {
   POLY_CHECK( position.size() == 3 );
   BGPoint<RealType,3> point(position[0], position[1], position[2]);
   return point;
}

//------------------------------------------------------------------------------
// Make a Boost.Geometry polygon from a concatenated vector of (x,y) points
//------------------------------------------------------------------------------
template <typename RealType>
BGPolygon<RealType, 2>
makeBGPolygon( std::vector<RealType>& points ) {
  BGPolygon<RealType, 2> polygon;
   for (unsigned i = 0; i < points.size()/2; ++i) {
     boost::geometry::append( polygon, BGPoint<RealType,2>(points[2*i],points[2*i+1]) );
   }
   boost::geometry::append( polygon, BGPoint<RealType,2>(points[0],points[1]) );
   return polygon;
}

//------------------------------------------------------------------------------
// Make a Boost.Geometry polygon from a PLC and its point list
//------------------------------------------------------------------------------
template <typename RealType>
BGPolygon<RealType, 2>
makeBGPolygon( PLC<2>& PLC, std::vector<RealType>& PLCpoints ) {
   unsigned i,j;
   BGPolygon<RealType, 2> polygon;
   // Walk the facets and add the first node
   for (j = 0; j != PLC.facets.size(); ++j){
      POLY_CHECK( PLC.facets[j].size() == 2 );
      i = PLC.facets[j][0];
      boost::geometry::append( polygon, BGPoint<RealType,2>(PLCpoints[2*i],PLCpoints[2*i+1]) );
   }
   i = PLC.facets[0][0];
   boost::geometry::append( polygon, BGPoint<RealType,2>(PLCpoints[2*i],PLCpoints[2*i+1]) ); //Close the polygon
   
   // Walk the facets composing each hole and add the first node
   const unsigned nHoles = PLC.holes.size();
   if (nHoles > 0) {
      auto& holes = polygon.inners();
      holes.resize(nHoles);
      for (unsigned ihole = 0; ihole != nHoles; ++ihole) {
         for (j = 0; j != PLC.holes[ihole].size(); ++j){
            POLY_CHECK( PLC.holes[ihole][j].size() == 2 );
            i = PLC.holes[ihole][j][0];
            boost::geometry::append( holes[ihole], BGPoint<RealType,2>( PLCpoints[2*i], PLCpoints[2*i+1] ) );
         }
         i = PLC.holes[ihole][0][0];
         boost::geometry::append( holes[ihole], BGPoint<RealType,2>( PLCpoints[2*i], PLCpoints[2*i+1] ) );  //Close the polygon
      }
   }
   return polygon;
}

}

#endif
