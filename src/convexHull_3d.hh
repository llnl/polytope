//----------------------------------------------------------------------------//
// 3D implementation of the convex hull algorithm.
// My own (possibly loose) interpretation of the incremental hull algorithm
// in 3D (JMO).
//----------------------------------------------------------------------------//
#ifndef __polytope_convexHull_3d__
#define __polytope_convexHull_3d__

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>

#include "PLC.hh"
#include "polytope_internal.hh"
#include "polytope_geometric_utilities.hh"
#include "Point.hh"

namespace polytope {

namespace {

//------------------------------------------------------------------------------
// Distance from a point to a line.
// The line is specified by a point on the line (line0) and a unit direction
// vector (lineDirection).
//------------------------------------------------------------------------------
template<typename RealType>
RealType
pointLineDistance(const Point3<RealType>& p,
                  const Point3<RealType>& line0,
                  const Point3<RealType>& lineDirection) {
  const Point3<RealType> delta = p - line0;
  const RealType a = geometry::dot<3, RealType>(&delta.x, &lineDirection.x);
  const RealType c2 = geometry::dot<3, RealType>(&delta.x, &delta.x);
  return sqrt(std::abs(c2 - a*a));
}

//------------------------------------------------------------------------------
// Signed distance from a point to a plane
// The plane is specified by a point on the plane (line0) and a unit normal to
// the plane (planeNormal).
// The return value is signed based on if the point is above (> 0) or below (< 0)
// the plane based on the normal direction.
//------------------------------------------------------------------------------
template<typename RealType>
RealType
pointPlaneDistance(const Point3<RealType>& p,
                   const Point3<RealType>& plane0,
                   const Point3<RealType>& planeNormal) {
  const Point3<RealType> delta = p - plane0;
  return geometry::dot<3, RealType>(&delta.x, &planeNormal.x);
}

//------------------------------------------------------------------------------
// Comparison by (x,y, or z) coordinate.
//------------------------------------------------------------------------------
template<typename RealType>
bool compareMinX(const Point3<RealType>& a, const Point3<RealType>& b) {
  return a.x < b.x;
}

template<typename RealType>
bool compareMinY(const Point3<RealType>& a, const Point3<RealType>& b) {
  return a.y < b.y;
}

template<typename RealType>
bool compareMinZ(const Point3<RealType>& a, const Point3<RealType>& b) {
  return a.z < b.z;
}

//------------------------------------------------------------------------------
// Compare points by distance from a line.
//------------------------------------------------------------------------------
template<typename RealType>
struct CompareLineDistance {
  Point3<RealType> point, direction;
  CompareLineDistance(const Point3<RealType>& p,
                      const Point3<RealType>& dir): point(p), direction(dir) {}
  bool operator()(const Point3<RealType>& a, const Point3<RealType>& b) {
    const RealType adist = pointLineDistance(a, point, direction),
                   bdist = pointLineDistance(b, point, direction);
    return adist < bdist;
  }
};

//------------------------------------------------------------------------------
// Compare points by distance from a plane.
//------------------------------------------------------------------------------
template<typename RealType>
struct CompareAbsPlaneDistance {
  Point3<RealType> point, normal;
  CompareAbsPlaneDistance(const Point3<RealType>& p,
                          const Point3<RealType>& norm): point(p), normal(norm) {}
  bool operator()(const Point3<RealType>& a, const Point3<RealType>& b) {
    const RealType adist = std::abs(pointPlaneDistance(a, point, normal)),
                   bdist = std::abs(pointPlaneDistance(b, point, normal));
    return adist < bdist;
  }
};

// //------------------------------------------------------------------------------
// // Borrow the Point3 type as a tuple to create 3 node facets hashes.
// //------------------------------------------------------------------------------
// Point3<unsigned>
// hashFacet(const unsigned i, const unsigned j, const unsigned k) {
//   typedef Point3<unsigned> Tuple3;
//   POLY_ASSERT(i != j and i != k and j != k);
//   if (i < j and i < k) {
//     if (j < k) {
//       return Tuple3(i, j, k);
//     } else {
//       return Tuple3(i, k, j);
//     }
//   } else if (j < i and j < k) {
//     if (i < k) {
//       return Tuple3(j, i, k);
//     } else {
//       return Tuple3(j, k, i);
//     }
//   } else {
//     if (i < j) {
//       return Tuple3(k, i, j);
//     } else {
//       return Tuple3(k, j, i);
//     }
//   }
// }

//------------------------------------------------------------------------------
// Hold one of our triangular facets.
//------------------------------------------------------------------------------
template<typename RealType>
struct TriangleFacet {
  typedef Point3<RealType> PointType;
  typedef Point3<unsigned> HashType;

  unsigned inode, jnode, knode;          // Indices of the vertices of the triangle.
  PointType normal;                      // Normal to the surface.

  TriangleFacet(): inode(0), jnode(0), knode(0), normal() {}
  TriangleFacet(const unsigned i, 
                const unsigned j,
                const unsigned k,
                const std::vector<PointType>& points): inode(i), 
                                                       jnode(j),
                                                       knode(k),
                                                       normal() {
    this->computeNormal(points);
  }

  void computeNormal(const std::vector<PointType>& points) {
    const PointType ij = points[jnode] - points[inode];
    const PointType ik = points[knode] - points[inode];
    geometry::cross<3, RealType>(&ij.x, &ik.x, &normal.x);
    geometry::unitVector<3, RealType>(&normal.x);
  }

  int compare(const PointType& point, 
              const std::vector<PointType>& points) const { 
    const PointType ip = point - points[inode];
    const RealType test = geometry::dot<3, RealType>(&ip.x, &normal.x);
    return (std::abs(test) < 1.0e-8 ? 0 :
            test < 0 ? -1 :
            1);
  }

  void flip() { std::swap(inode, jnode); normal *= -1; }
  // bool operator==(const TriangleFacet& rhs) const { return hashFacet(inode, jnode, knode) == hashFacet(rhs.inode, rhs.jnode, rhs.knode); }
  // bool operator<(const TriangleFacet& rhs) const { return hashFacet(inode, jnode, knode) < hashFacet(rhs.inode, rhs.jnode, rhs.knode); }

  // output operator
  friend std::ostream& operator<<(std::ostream& os, const TriangleFacet& facet) {
    os << "TriangleFacet(" << facet.inode << " " << facet.jnode << " " << facet.knode 
       << "), normal = " << facet.normal;
    return os;
  }
};

//------------------------------------------------------------------------------
// Reduce the given set of point indices to those exterior to the facets.
//------------------------------------------------------------------------------
template<typename RealType>
void 
exteriorPoints(std::vector<unsigned>& indices,
               const std::vector<TriangleFacet<RealType> >& facets, 
               const std::vector<Point3<RealType> >& points) {
  typedef std::vector<TriangleFacet<RealType> > FacetSet;

  std::vector<unsigned> remainingPoints;
  for (std::vector<unsigned>::const_iterator iitr = indices.begin();
       iitr != indices.end();
       ++iitr) {
    typename FacetSet::const_iterator fitr = facets.begin();
    while (fitr != facets.end() and fitr->compare(points[*iitr], points) != 1) ++fitr;
    if (fitr != facets.end()) remainingPoints.push_back(*iitr);
  }
  indices = remainingPoints;
}

//------------------------------------------------------------------------------
// Reduce the given set of point indices to those exterior to the facets.
//------------------------------------------------------------------------------
template<typename RealType>
unsigned highestPoint(std::vector<unsigned>& indices,
                      const std::vector<TriangleFacet<RealType> >& facets, 
                      const std::vector<Point3<RealType> >& points) {
  typedef std::vector<TriangleFacet<RealType> > FacetSet;
  int result = -1;
  RealType maxAltitude = 0.0;
  for (typename std::vector<unsigned>::const_iterator itr = indices.begin();
       itr != indices.end();
       ++itr) {
    for (typename FacetSet::const_iterator fitr = facets.begin();
         fitr != facets.end();
         ++fitr) {
      const RealType altitude = pointPlaneDistance(points[*itr], points[fitr->inode], fitr->normal);
      if (altitude > maxAltitude) {
        maxAltitude = altitude;
        result = *itr;
      }
    }
  }
  POLY_ASSERT(result >= 0);
  return unsigned(result);
}

//------------------------------------------------------------------------------
// Cull out the facets visible from the given point (points[apex]).
// Also computes the remaining horizon edges.
//------------------------------------------------------------------------------
template<typename RealType>
void cullVisibleFacets(std::vector<std::pair<int, int> >& horizonEdges,
                       std::vector<TriangleFacet<RealType> >& facets, 
                       const std::vector<Point3<RealType> >& points,
                       const unsigned apex) {

  typedef std::pair<int, int> EdgeHash;
  typedef std::vector<TriangleFacet<RealType> > FacetSet;

  FacetSet newFacets;
  internal::CounterMap<EdgeHash> edgeCount;

  // Walk the facets.
  for (typename FacetSet::const_iterator fitr = facets.begin();
       fitr != facets.end();
       ++fitr) {

    // Is this facet visible from the apex point?
    if (fitr->compare(points[apex], points) == 1) {

      // Yep, it needs to be removed.  Increment the use count for each
      // of it's edges.
      ++edgeCount[internal::hashEdge(fitr->inode, fitr->jnode)];
      ++edgeCount[internal::hashEdge(fitr->jnode, fitr->knode)];
      ++edgeCount[internal::hashEdge(fitr->knode, fitr->inode)];

    } else {

      // Nope, copy it to the new set.
      newFacets.push_back(*fitr);

    }
  }
  facets = newFacets;

  // The horizon edges are the ones with a single use count.  All others should
  // have been hit twice.
  horizonEdges = std::vector<EdgeHash>();
  for (internal::CounterMap<EdgeHash>::const_iterator eitr = edgeCount.begin();
       eitr != edgeCount.end();
       ++eitr) {
    POLY_ASSERT(eitr->second == 1 or eitr->second == 2);
    if (eitr->second == 1) horizonEdges.push_back(eitr->first);
  }
}

} // anonymous namespace

//------------------------------------------------------------------------------
// The 3D convex hull method.
//------------------------------------------------------------------------------
template<typename RealType>
PLC<3>
convexHull_3d(const std::vector<RealType>& points,
              const RealType* low,
              const RealType& dx) {

  typedef Point3<RealType> PointType;
  typedef std::vector<TriangleFacet<RealType> > FacetSet;
  typedef std::pair<int, int> EdgeHash;

  // Pre-conditions.
  POLY_ASSERT(points.size() % 3 == 0);

  // Compute the bounding box.
  RealType xmin[3], xmax[3];
  geometry::computeBoundingBox<3, RealType>(points, true, xmin, xmax);

  // Reduce to the unique set of input points.
  std::vector<RealType> upoints;
  std::vector<unsigned> pointMap;
  geometry::uniquePoints<3, RealType>(points, xmin, xmax, dx, upoints, pointMap);
  const unsigned nunique = upoints.size()/3;
  POLY_ASSERT(nunique >= 4);
  std::vector<PointType> ps(nunique);
  for (unsigned i = 0; i != nunique; ++i) {
    ps[i].x = upoints[3*i];
    ps[i].y = upoints[3*i+1];
    ps[i].z = upoints[3*i+2];
    ps[i].index = i;
  }

  // Find the points on the min/max x, y, & z.  These must be in the hull.
  const unsigned ixmin = std::min_element(ps.begin(), ps.end(), compareMinX<RealType>)->index,
                 ixmax = std::max_element(ps.begin(), ps.end(), compareMinX<RealType>)->index,
                 iymin = std::min_element(ps.begin(), ps.end(), compareMinY<RealType>)->index,
                 iymax = std::max_element(ps.begin(), ps.end(), compareMinY<RealType>)->index,
                 izmin = std::min_element(ps.begin(), ps.end(), compareMinZ<RealType>)->index,
                 izmax = std::max_element(ps.begin(), ps.end(), compareMinZ<RealType>)->index;
  POLY_ASSERT(ixmin != ixmax and
              iymin != iymax and
              izmin != izmax);

  // Choose the two points furthest apart.
  const PointType boxx = ps[ixmax] - ps[ixmin], boxy = ps[iymax] - ps[iymin], boxz = ps[izmax] - ps[izmin];
  const RealType boxx2 = geometry::dot<3, RealType>(&boxx.x, &boxx.x), 
                 boxy2 = geometry::dot<3, RealType>(&boxy.x, &boxy.x),
                 boxz2 = geometry::dot<3, RealType>(&boxz.x, &boxz.x);
  POLY_ASSERT(boxx2 > 0.0 and boxy2 > 0.0 and boxz2 > 0.0);
  TriangleFacet<RealType> startingFacet;
  if (boxx2 >= std::max(boxy2, boxz2)) {
    startingFacet.inode = ixmin;
    startingFacet.jnode = ixmax;
  } else if (boxy2 >= boxz2) {
    startingFacet.inode = iymin;
    startingFacet.jnode = iymax;
  } else {
    POLY_ASSERT(boxz2 >= std::max(boxx2, boxy2));
    startingFacet.inode = izmin;
    startingFacet.jnode = izmax;
  }

  // Select a third point as the most distant from the line defined by the two we just picked.
  // This one should also be on the hull.
  PointType lineDirection = ps[startingFacet.jnode] - ps[startingFacet.inode];
  geometry::unitVector<3, RealType>(&lineDirection.x);
  startingFacet.knode = std::max_element(ps.begin(), ps.end(), CompareLineDistance<RealType>(ps[startingFacet.inode], lineDirection))->index;
  startingFacet.computeNormal(ps);

  // We have the triangular base, so pick the furthest point from this base as the apex of our
  // starting tetrahedron.  This point should also be on the hull.
  const unsigned apex = std::max_element(ps.begin(), ps.end(), CompareAbsPlaneDistance<RealType>(ps[startingFacet.inode], startingFacet.normal))->index;

  // If the apex point is below our starting facet, flip the facet.
  if (startingFacet.compare(ps[apex], ps) == 1) startingFacet.flip();

  // Intialize our starting tetrahedron.
  FacetSet workingFacets;
  workingFacets.push_back(startingFacet);
  workingFacets.push_back(TriangleFacet<RealType>(startingFacet.jnode, startingFacet.inode, apex, ps));
  workingFacets.push_back(TriangleFacet<RealType>(startingFacet.inode, startingFacet.knode, apex ,ps));
  workingFacets.push_back(TriangleFacet<RealType>(startingFacet.knode, startingFacet.jnode, apex, ps));
  // std::cerr << "Starting facets: " << std::endl
  //           << "    " << workingFacets[0] << std::endl
  //           << "    " << workingFacets[1] << std::endl
  //           << "    " << workingFacets[2] << std::endl
  //           << "    " << workingFacets[3] << std::endl;

  // Find the set of points exterior to our starting tetrahedron.
  std::vector<unsigned> remainingPoints(ps.size());
  for (unsigned i = 0; i != ps.size(); ++i) remainingPoints[i] = i;
  exteriorPoints(remainingPoints, workingFacets, ps);

  // Iterate until all points are interior to the facets.
  while (!remainingPoints.empty()) {

    // Find the remaining point furthest from the hull.
    const unsigned apex = highestPoint(remainingPoints, workingFacets, ps);

    // Remove the visible facets from apex, and compute the horizon edges.
    std::vector<EdgeHash> horizonEdges;
    cullVisibleFacets(horizonEdges, workingFacets, ps, apex);
    POLY_ASSERT(horizonEdges.size() >= 3);

    // Build the new facets from the apex point and horizon edges.
    // const PointType testPoint = (ps[workingFacets[0].inode] + ps[workingFacets[0].jnode] + ps[workingFacets[0].knode])/3.0;
    PointType testPoint;
    for (typename std::vector<EdgeHash>::const_iterator eitr = horizonEdges.begin();
         eitr != horizonEdges.end();
         ++eitr) {
      testPoint += ps[eitr->first];
      testPoint += ps[eitr->second];
    }
    testPoint /= 2*horizonEdges.size();
    for (typename std::vector<EdgeHash>::const_iterator eitr = horizonEdges.begin();
         eitr != horizonEdges.end();
         ++eitr) {
      workingFacets.push_back(TriangleFacet<RealType>(eitr->first, eitr->second, apex, ps));
      if (workingFacets.back().compare(testPoint, ps) == 1) workingFacets.back().flip();
    }

    // Reduce to the new set of exterior points.
    const unsigned oldsize = remainingPoints.size();
    POLY_CONTRACT_VAR(oldsize);
    exteriorPoints(remainingPoints, workingFacets, ps);
    POLY_ASSERT(remainingPoints.size() < oldsize);

    // std::cerr << "New pass at facets: " << std::endl;
    // for (unsigned i = 0; i != workingFacets.size(); ++i) {
    //   std::cerr << " --> " << workingFacets[i] << std::endl;
    // }
    // std::cerr << "Remaining points : ";
    // std::copy(remainingPoints.begin(), remainingPoints.end(), std::ostream_iterator<unsigned>(std::cerr, " "));
    // std::cerr << std::endl;
    // std::cerr << "  --> " << remainingPoints.size() << std::endl;
  }

  // Read out the data to the PLC and we're done.
  PLC<3> plc;
  for (typename FacetSet::iterator facetItr = workingFacets.begin();
       facetItr != workingFacets.end();
       ++facetItr) {
    plc.facets.push_back(std::vector<int>());
    plc.facets.back().push_back(pointMap[facetItr->inode]);
    plc.facets.back().push_back(pointMap[facetItr->jnode]);
    plc.facets.back().push_back(pointMap[facetItr->knode]);
  }
  return plc;
}

}

#endif
