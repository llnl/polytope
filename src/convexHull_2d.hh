//----------------------------------------------------------------------------//
// 2D implementation of the convex hull algorithm.
// Based on an example at http://www.algorithmist.com/index.php/Monotone_Chain_Convex_Hull.cpp
//----------------------------------------------------------------------------//
#ifndef __polytope_convexHull_2d__
#define __polytope_convexHull_2d__

#include <iostream>
#include <iterator>
#include <algorithm>
#include <set>
#include <stdint.h>

#include "PLC.hh"
#include "polytope_internal.hh"
#include "polytope_geometric_utilities.hh"
#include "Point.hh"

namespace polytope {

namespace { // We hide internal functions in an anonymous namespace.

//------------------------------------------------------------------------------
// A fuzzy comparison operator for our quantized Point2 type.
//------------------------------------------------------------------------------
template<typename UintType>
struct FuzzyPoint2LessThan {
  UintType fuzz;
  FuzzyPoint2LessThan(const UintType ifuzz = 1): fuzz(ifuzz) {}
  bool operator()(const Point2<UintType>& p1, const Point2<UintType>& p2) const {
    return (p1.x + fuzz < p2.x        ? true :
            p1.x        > p2.x + fuzz ? false :
            p1.y + fuzz < p2.y        ? true :
            p1.y        > p2.y + fuzz ? false :
            false);
  }
  bool operator()(const std::pair<Point2<UintType>, unsigned>& p1,
                  const std::pair<Point2<UintType>, unsigned>& p2) const {
    return operator()(p1.first, p2.first);
  }
};

//------------------------------------------------------------------------------
// sign of the Z coordinate of cross product : (p2 - p1)x(p3 - p1).
//------------------------------------------------------------------------------
template<typename RealType>
int zcross_sign(const Point2<RealType>& p1, const Point2<RealType>& p2, const Point2<RealType>& p3) {
//   double scale = 1.0/max(RealType(1), max(p1.x, max(p1.y, max(p2.x, max(p2.y, max(p3.x, p3.y))))));
  const double ztest = 
    (double(p2.x) - double(p1.x))*(double(p3.y) - double(p1.y)) -
    (double(p2.y) - double(p1.y))*(double(p3.x) - double(p1.x));
  return (ztest < 0.0 ? -1 :
          ztest > 0.0 ?  1 :
                         0);
  // return (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x);
}

//------------------------------------------------------------------------------
// Comparator to compare std::pair's by their first element.
//------------------------------------------------------------------------------
template<typename T1, typename T2>
struct ComparePairByFirstElement {
  bool operator()(const std::pair<T1, T2>& lhs, const std::pair<T1, T2>& rhs) const {
    return lhs.first < rhs.first;
  }
};

} // end anonymous namespace

//------------------------------------------------------------------------------
// The method itself.
//
// NOTE: The convex hull can be of dimension smaller than 2D. Lower
//       dimensionality is stored in the structure of the PLC facets
//             1D - Collinear points - A single length-2 facet with indices 
//                                     pointing to the smallest and largest
//                                     points in the sorted point hash
//             0D - Single point     - A single length-2 facet with the same
//                                     index (0) in both positions 
//------------------------------------------------------------------------------
template<typename RealType>
PLC<2, RealType>
convexHull_2d(const std::vector<RealType>& points,
              const RealType* low,
              const RealType& dx) {
  typedef KeyTraits::Key CoordHash;
  typedef Point2<CoordHash> PointHash;
  // typedef polytope::DimensionTraits<2, RealType>::CoordHash CoordHash;
  // typedef polytope::DimensionTraits<2, RealType>::IntPoint PointHash;

  POLY_ASSERT(!points.empty());
  POLY_ASSERT(points.size() % 2 == 0);
  const unsigned n = points.size() / 2;
  PLC<2, RealType> plc;
  unsigned i, j, k, t;
  
  // If there's only one or two points, we're done: that's the whole hull
  if (n == 1 or n == 2) {
    plc.facets.resize(1, std::vector<int>(2));
    plc.facets[0][0] = 0;
    plc.facets[0][1] = (n == 1) ? 0 : 1;
    return plc;
  }
  
  // Start by finding a point distinct from point 0.
  j = 1;
  while (j != n and geometry::distance<2, RealType>(&points[0], &points[2*j]) < dx) ++j;
  if (j == n - 1) {
    // There are only 2 distinct positions!
    plc.facets.resize(1, std::vector<int>(2));
    plc.facets[0][0] = 0;
    plc.facets[0][1] = j;
    return plc;
  } else if (j == n) {
    // Good god, there are no distinct points!
    plc.facets.resize(1, std::vector<int>(2));
    plc.facets[0][0] = 0;
    plc.facets[0][1] = 0;
    return plc;
  }

  // Check if the input points are collinear.
  bool collinear = true;
  POLY_ASSERT(n > 2);
  i = 2;
  while (collinear and i != n) {
    collinear = geometry::collinear<2,RealType>(&points[0], &points[2*j], &points[2*i], dx);
    ++i;
  }
  
  // Hash the input points and sort them by x coordinate, remembering their original indices
  // in the input set.  We also ensure that only unique (using a fuzzy comparison) points
  // are inserted here, since duplicates mess up the hull calculation.
  const RealType& xmin = low[0];
  const RealType& ymin = low[1];
  std::set<std::pair<PointHash, unsigned>, FuzzyPoint2LessThan<CoordHash> > uniquePoints;
  for (i = 0; i != n; ++i) {
    uniquePoints.insert(std::make_pair(PointHash(CoordHash((points[2*i]     - xmin)/dx + 0.5),
                                                 CoordHash((points[2*i + 1] - ymin)/dx + 0.5)),
                                       i));
  }
  std::vector<std::pair<PointHash, unsigned> > sortedPoints(uniquePoints.begin(), uniquePoints.end());
  std::sort(sortedPoints.begin(), sortedPoints.end());

  // If the points are collinear, we can save a lot of work
  if (collinear) {
    plc.facets.resize(1, std::vector<int>(2));
    plc.facets[0][0] = sortedPoints.front().second;
    plc.facets[0][1] = sortedPoints.back().second;
  }
  else {
    // Prepare the result.
    const unsigned nunique = sortedPoints.size();
    std::vector<int> result(2*nunique);
    
    // Build the lower hull.
    for (i = 0, k = 0; i < nunique; i++) {
      while (k >= 2 and
             zcross_sign(sortedPoints[result[k - 2]].first, sortedPoints[result[k - 1]].first, sortedPoints[i].first) <= 0) k--;
      result[k++] = i;
    }
    
    // Build the upper hull.
    for (i = nunique - 2, t = k + 1; i >= 0; i--) {
      while (k >= t and
             zcross_sign(sortedPoints[result[k - 2]].first, sortedPoints[result[k - 1]].first, sortedPoints[i].first) <= 0) k--;
      result[k++] = i;
    }
    // if (!(k >= 4)) {
    //   std::cerr << "Blago!  " << n << " " << nunique << " " << k << std::endl;
    //   std::cerr << "Unique:" << std::endl;
    //   for (unsigned i = 0; i != nunique; ++i) std::cerr << "  --> " << sortedPoints[i].first << std::endl;
    //   std::cerr << "Input:" << std::endl;
    //   for (unsigned i = 0; i != n; ++i) std::cerr << "  --> " << points[2*i] << " " << points[2*i+1] << std::endl;
    // }
    POLY_ASSERT(k >= 4);
    POLY_ASSERT(result.front() == result.back());
    
    // Translate our sorted information to a PLC based on the input point ordering and we're done.
    for (i = 0; i != k - 1; ++i) {
      j = (i + 1) % k;
      plc.facets.push_back(std::vector<int>());
      plc.facets.back().push_back(sortedPoints[result[i]].second);
      plc.facets.back().push_back(sortedPoints[result[j]].second);
    }
    POLY_ASSERT(plc.facets.size() == k - 1);
  }
  return plc;
}

}

#endif
