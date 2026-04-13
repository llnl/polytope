#ifndef POLYTOPE_PLC_BOOST_2D_HH
#define POLYTOPE_PLC_BOOST_2D_HH

#ifdef POLYTOPE_ENABLE_BOOST

#include <cstddef>
#include <vector>

#include "ReducedPLC.hh"
#include "Point.hh"

// Computational geometetry operations use Boost.Geometry library.
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

namespace polytope {
namespace BG {

//------------------------------------------------------------------------------
// Public interface methods (forward declaration).
//------------------------------------------------------------------------------
  template<typename RealType> std::vector<ReducedPLC<2, RealType> > boost_union    (const ReducedPLC<2, RealType>& a, const ReducedPLC<2, RealType>& b);
  template<typename RealType> std::vector<ReducedPLC<2, RealType> > boost_intersect(const ReducedPLC<2, RealType>& a, const ReducedPLC<2, RealType>& b);
  template<typename RealType> std::vector<ReducedPLC<2, RealType> > boost_subtract (const ReducedPLC<2, RealType>& a, const ReducedPLC<2, RealType>& b);
  template<typename RealType> ReducedPLC<2, RealType> boost_clip(const ReducedPLC<2, RealType>& a,
								 const ReducedPLC<2, RealType>& b,
								 const Point2<RealType>& p,
								 std::vector<ReducedPLC<2, RealType> >& orphans);

  template<typename RealType> bool boost_within(const Point2<RealType>& p, const ReducedPLC<2, RealType>& a);
  template<typename RealType> bool boost_intersects(const ReducedPLC<2, RealType>& a, const ReducedPLC<2, RealType>& b);
  template<typename RealType> bool boost_intersects(const ReducedPLC<2, RealType>& a);

   template<typename RealType> std::vector<ReducedPLC<2, RealType> > boost_unionReduce(const std::vector<ReducedPLC<2, RealType> >& a);


//------------------------------------------------------------------------------
// Internal methods called by boost
//------------------------------------------------------------------------------
namespace PLC_boost_internal {

// Convert a ReducedPLC to a boost::geometry polygon
template<typename RealType>
boost::geometry::model::polygon<Point2<RealType>, false>
ReducedPLCtoPolygon(const ReducedPLC<2, RealType>& model) {
  typedef Point2<RealType> PointType;
  typedef boost::geometry::model::polygon<PointType, false> PolygonType;
  const unsigned nfacets = model.facets.size();
  const unsigned nholes  = model.holes.size();
  PolygonType result;
  int ii = model.facets[0][0];
  boost::geometry::append(result, 
			  PointType(model.points[2*ii], model.points[2*ii+1]));
  for (unsigned j = 0; j != nfacets; ++j) {
    POLY_ASSERT(model.facets[j].size() == 2);
    ii = model.facets[j][1];
    boost::geometry::append(result, 
			    PointType(model.points[2*ii], model.points[2*ii+1]));
  }
  POLY_ASSERT(result.outer().size()  == nfacets + 1);
  POLY_ASSERT(result.outer().front() == result.outer().back());
  
  typename PolygonType::inner_container_type& holes = result.inners();
  holes.resize(nholes);
  for (unsigned j = 0; j != nholes; ++j) {
    const unsigned nhfacets = model.holes[j].size();
    ii = model.holes[j][0][0];
    boost::geometry::append(holes[j], 
			    PointType(model.points[2*ii], model.points[2*ii+1]));
    for (unsigned k = 0; k != nhfacets; ++k) {
      POLY_ASSERT(model.holes[j][k].size() == 2);
      ii = model.holes[j][k][1];
      boost::geometry::append(holes[j], 
			      PointType(model.points[2*ii], model.points[2*ii+1]));
    }
    POLY_ASSERT(holes[j].size()  == nhfacets + 1);
    POLY_ASSERT(holes[j].front() == holes[j].back());
  }
  POLY_ASSERT(result.inners().size() == nholes);
  boost::geometry::correct(result);
  return result;
}


// Convert a boost::geometry polygon to a ReducedPLC
template<typename RealType>
ReducedPLC<2, RealType>
ReducedPLCfromPolygon(const boost::geometry::model::polygon<Point2<RealType>, false>& model) {
  typedef Point2<RealType> PointType;
  typedef boost::geometry::model::polygon<PointType, false> PolygonType;
  typedef boost::geometry::model::ring<PointType, false> RingType;
  const unsigned nfacets = model.outer().size() - 1;
  const unsigned nholes  = model.inners().size();
  POLY_ASSERT(nfacets >= 3);
  ReducedPLC<2, RealType> result;
  int offset = 0;
  int ii = 0;
  // Add the facets for the outer boundary
  result.facets.resize(nfacets, std::vector<int>(2));
  for (typename RingType::const_iterator itr = model.outer().begin();
       itr != model.outer().end() - 1;
       ++itr, ii++) {
    result.points.push_back((*itr).x);
    result.points.push_back((*itr).y);
    result.facets[ii][0] = ii;
    result.facets[ii][1] = (ii+1) % nfacets;
  }
  POLY_ASSERT2(ii == nfacets, ii << " " << nfacets);
  offset += ii;
  // Now add in any holes
  result.holes.resize(nholes);
//   const typename PolygonType::inner_container_type& allHoles = model.inners();
  for (unsigned ihole = 0; ihole != nholes; ++ihole) {
    const typename PolygonType::ring_type& hole = model.inners()[ihole];
    const unsigned nhfacets = hole.size() - 1;
    result.holes[ihole].resize(nhfacets, std::vector<int>(2));
    ii = 0;
    for (typename PolygonType::ring_type::const_iterator itr = hole.begin();
//     for (typename PolygonType::inner_container_type::const_iterator itr = hole.begin();
	 itr != hole.end() - 1;
	 ++itr, ii++) {
      result.points.push_back((*itr).x);
      result.points.push_back((*itr).y);
      result.holes[ihole][ii][0] = offset + ii;
      result.holes[ihole][ii][0] = offset + ((ii+1) % nhfacets);
    }
    POLY_ASSERT(ii == nhfacets-1);
    offset += ii;
  }
  return result;
}


// // Determine the single polygon in a collection that contains the given point.
// // Return the index of the polygon in the collection
// template<typename RealType>
// int polygonContainingPoint(const std::vector<boost::geometry::model::polygon<Point2<RealType>, false> > polygons,
// 			   const Point2<RealType> point) {
//   if (polygons.size() == 1) {
//     POLY_ASSERT(boost::geometry::covered_by(point, intersections[0]));
//     return 0;
//   } else {
//     bool inside = false;
//     unsigned j = 0;
//     while (j != polygons.size() and not inside) {
//       inside = boost::geometry::covered_by(point, polygons[j]);
//       if (not inside) {
// 	for (typename RingType::const_iterator itr = polygons[j].outer().begin();
// 	     itr != polygons[j].outer().end() - 1;
// 	     ++itr)  inside += (point == *itr);
//       }
//       j++;
//     }
//     POLY_ASSERT2(j != polygons.size(),
// 		 "None of the input polygons contain point " << point);
//     return j-1;
//   }
// }

} // end namespace PLC_boost_internal



//------------------------------------------------------------------------------
// Public interface implementations
//------------------------------------------------------------------------------

// Union
template<typename RealType>
std::vector<ReducedPLC<2, RealType> > boost_union(const ReducedPLC<2, RealType>& a,
						  const ReducedPLC<2, RealType>& b) {
  using namespace PLC_boost_internal;
  typedef boost::geometry::model::polygon<Point2<RealType>, false> PolygonType;
  std::vector<PolygonType> bgresult;
  boost::geometry::union_(ReducedPLCtoPolygon(a), ReducedPLCtoPolygon(b), bgresult);
  POLY_ASSERT(not bgresult.empty());
  std::vector<ReducedPLC<2, RealType> > result(bgresult.size());
  for (std::size_t i = 0; i != bgresult.size(); ++i) {
    POLY_ASSERT(bgresult[i].outer().front() == bgresult[i].outer().back());
    POLY_ASSERT(not boost::geometry::intersects(bgresult[i]));
    result[i] = ReducedPLCfromPolygon(bgresult[i]);
  }
  return result;
}


// Intersection
template<typename RealType>
std::vector<ReducedPLC<2, RealType> > boost_intersect(const ReducedPLC<2, RealType>& a,
						      const ReducedPLC<2, RealType>& b) {
  using namespace PLC_boost_internal;
  typedef boost::geometry::model::polygon<Point2<RealType>, false> PolygonType;
  std::vector<PolygonType> bgresult;
  boost::geometry::intersection(ReducedPLCtoPolygon(a), ReducedPLCtoPolygon(b), bgresult);
  POLY_ASSERT(not bgresult.empty());
  std::vector<ReducedPLC<2, RealType> > result(bgresult.size());
  for (std::size_t i = 0; i != bgresult.size(); ++i) {
    POLY_ASSERT(bgresult[i].outer().front() == bgresult[i].outer().back());
    POLY_ASSERT(not boost::geometry::intersects(bgresult[i]));
    result[i] = ReducedPLCfromPolygon(bgresult[i]);
  }
  return result;
}


// Subtraction
template<typename RealType>
std::vector<ReducedPLC<2, RealType> > boost_subtract(const ReducedPLC<2, RealType>& a,
						     const ReducedPLC<2, RealType>& b) {
  using namespace PLC_boost_internal;
  typedef boost::geometry::model::polygon<Point2<RealType>, false> PolygonType;
  std::vector<PolygonType> bgresult;
  boost::geometry::difference(ReducedPLCtoPolygon(a), ReducedPLCtoPolygon(b), bgresult);
  POLY_ASSERT(not bgresult.empty());
  std::vector<ReducedPLC<2, RealType> > result(bgresult.size());
  for (std::size_t i = 0; i != bgresult.size(); ++i) {
    POLY_ASSERT(bgresult[i].outer().front() == bgresult[i].outer().back());
    POLY_ASSERT(not boost::geometry::intersects(bgresult[i]));
    result[i] = ReducedPLCfromPolygon(bgresult[i]);
  }
  return result;
}


// Intersection, returning the single ReducedPLC containing point p
template<typename RealType>
ReducedPLC<2, RealType> boost_clip(const ReducedPLC<2, RealType>& a,
				   const ReducedPLC<2, RealType>& b,
				   const Point2<RealType>& p,
				   std::vector<ReducedPLC<2, RealType> >& orphans) {
  using namespace PLC_boost_internal;
  typedef boost::geometry::model::ring<Point2<RealType>, false> RingType;
  typedef boost::geometry::model::polygon<Point2<RealType>, false> PolygonType;
  std::vector<PolygonType> intersections;
  boost::geometry::intersection(ReducedPLCtoPolygon(a), ReducedPLCtoPolygon(b), intersections);
  POLY_ASSERT2(intersections.size() > 0,
	       p.x << " " << p.y << std::endl 
	       << boost::geometry::dsv(ReducedPLCtoPolygon(a)) << std::endl
	       << boost::geometry::dsv(ReducedPLCtoPolygon(b)) << std::endl);
  PolygonType result;
  
  // Only one geometry in the intersection
  if (intersections.size() == 1) {
    result = intersections[0];
  }
  // The one geometry that contains point p is the returned intersection
  else {
    bool inside, onBoundary;
    std::size_t k = intersections.size();
    for (std::size_t j = 0; j != intersections.size(); ++j) {
      inside = boost::geometry::covered_by(p, intersections[j]);
      if (inside) k = j;
      else {
	onBoundary = false;
	for (typename RingType::const_iterator itr = intersections[j].outer().begin();
	     itr != intersections[j].outer().end() - 1;
	     ++itr)  onBoundary += (p == *itr);
	if (onBoundary) k = j;
	else            orphans.push_back(ReducedPLCfromPolygon(intersections[j]));
      }
    }
    POLY_ASSERT(k < intersections.size());
    result = intersections[k];
  }
  boost::geometry::correct(result);
  boost::geometry::unique(result);
  // POLY_ASSERT2(not boost::geometry::intersects(result),
  //              "Self-intersecting result:\n" << boost::geometry::dsv(result));
  return ReducedPLCfromPolygon(result);
}


// Within -- Point in Polygon
template<typename RealType>
bool boost_within(const Point2<RealType>& p,
                  const ReducedPLC<2, RealType>& a) {
  using namespace PLC_boost_internal;
  return boost::geometry::within(p, ReducedPLCtoPolygon(a));
}

// Intersection -- Polygon and Polygon
template<typename RealType>
bool boost_intersects(const ReducedPLC<2, RealType>& a,
                      const ReducedPLC<2, RealType>& b) {
  using namespace PLC_boost_internal;
  return boost::geometry::intersects(ReducedPLCtoPolygon(a), ReducedPLCtoPolygon(b));
}

// Intersection -- Polygon self intersection
template<typename RealType>
bool boost_intersects(const ReducedPLC<2, RealType>& a) {
  using namespace PLC_boost_internal;
  return boost::geometry::intersects(ReducedPLCtoPolygon(a));
}

// UnionReduce -- Reduce a collection of PLCs via unions
template<typename RealType>
std::vector<ReducedPLC<2, RealType> > boost_unionReduce(const std::vector<ReducedPLC<2, RealType> >& a) {
  using namespace PLC_boost_internal;
  typedef boost::geometry::model::polygon<Point2<RealType>, false> PolygonType;
  typedef boost::geometry::model::multi_polygon<PolygonType>       MultiPolygonType;
  MultiPolygonType bgresult;
  for (std::size_t i = 0; i < a.size(); ++i) {
    MultiPolygonType tmp;
    boost::geometry::union_(bgresult, ReducedPLCtoPolygon(a[i]), tmp);
    boost::geometry::correct(tmp);
    bgresult = tmp;
  }
  POLY_ASSERT(not bgresult.empty());
  POLY_ASSERT(bgresult.size() <= a.size());
  std::vector<ReducedPLC<2, RealType> > result(bgresult.size());
  for (std::size_t i = 0; i < bgresult.size(); ++i) {
    POLY_ASSERT(bgresult[i].outer().front() == bgresult[i].outer().back());
    POLY_ASSERT(not boost::geometry::intersects(bgresult[i]));
    result[i] = ReducedPLCfromPolygon(bgresult[i]);
  }
  POLY_ASSERT(not result.empty());
  POLY_ASSERT(result.size() <= a.size());
  return result;
}  


} //end namespace BG
} //end namespace polytope


#endif  // end #ifdef POLYTOPE_ENABLE_BOOST

#endif
