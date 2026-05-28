#ifndef __Polytope_DimensionTraits__
#define __Polytope_DimensionTraits__
//------------------------------------------------------------------------------
// DimensionTraits
//
// Dimensional specializations to interface with Polytope's geometric objects
//------------------------------------------------------------------------------
#include <vector>

#include "KeyTraits.hh"
#include "convexHull_2d.hh"
#include "convexHull_3d.hh"
#include "QuantizedTessellation2d.hh"
#include "QuantizedTessellation3d.hh"
#include "ReducedPLC.hh"

namespace polytope {

// Base struct
template<int Dimension, typename RealType> struct DimensionTraits {};

// 2D specialization
template<typename RealType>
struct DimensionTraits<2, RealType> {
  typedef ReducedPLC<2, RealType> ConvexHull;
  typedef KeyTraits::Key CoordHash;
  typedef Point2<CoordHash> IntPoint;
  typedef Point2<RealType> RealPoint;
  typedef QuantizedTessellation2d<CoordHash, RealType> QuantizedTessellation;

  static ConvexHull convexHull(const std::vector<RealType>& points, 
                               const RealType* low,
                               const RealType& dx) { 
    return ConvexHull(convexHull_2d(points, low, dx), points);
  }
  static IntPoint constructPoint(const RealType* ri,
                                 const RealType* rlow,
                                 const RealType& dx,
                                 const size_t i) {
    return IntPoint(ri[0], ri[1], 
                    rlow[0], rlow[1], 
                    dx, i);
  }
  static RealPoint constructPoint(const RealType* ri) {
     return RealPoint(ri[0], ri[1]);
  }
  static IntPoint faceCentroid(const Tessellation<2, RealType>& mesh,
                               const unsigned iface,
                               const RealType* rlow,
                               const RealType& dx) {
    POLY_ASSERT(iface < mesh.faces.size());
    POLY_ASSERT(mesh.faces[iface].size() == 2);
    const unsigned n1 = mesh.faces[iface][0], n2 = mesh.faces[iface][1];
    POLY_ASSERT(n1 < mesh.nodes.size()/2);
    POLY_ASSERT(n2 < mesh.nodes.size()/2);
    const IntPoint pn1 = constructPoint(&mesh.nodes[2*n1], rlow, dx, n1);
    const IntPoint pn2 = constructPoint(&mesh.nodes[2*n2], rlow, dx, n2);
    const IntPoint result((pn1.x >> 1) + (pn2.x >> 1),
                          (pn1.y >> 1) + (pn2.y >> 1),
                          iface);
    return result;
  }
  static RealType maxLength(const RealType* low, const RealType* high) {
    return std::max(high[0] - low[0], high[1] - low[1]);
  }
  static std::vector<RealType> extractCoords(const std::vector<RealType>& allCoords,
                                             const std::vector<unsigned>& indices) {
    std::vector<RealType> result;
    result.reserve(2*indices.size());
    for (std::vector<unsigned>::const_iterator itr = indices.begin();
         itr != indices.end();
         ++itr) {
      const unsigned i = *itr;
      POLY_ASSERT(i < allCoords.size()/2);
      result.push_back(allCoords[2*i]);
      result.push_back(allCoords[2*i + 1]);
    }
    return result;
  }
  static int hullDimension(const ConvexHull hull) {
    if (hull.facets.size() == 1) {
      POLY_ASSERT(hull.facets[0].size() == 2);
      return (hull.facets[0][0] == hull.facets[0][1]) ? 0 : 1;
    } else return 2;
  }
};

// 3D specialization
template<typename RealType>
struct DimensionTraits<3, RealType> {
  typedef typename polytope::ReducedPLC<3, RealType> ConvexHull;
  typedef int64_t CoordHash;
  typedef polytope::Point3<CoordHash> IntPoint;
  typedef polytope::Point3<RealType> RealPoint;
  typedef polytope::QuantizedTessellation3d<CoordHash, RealType> QuantizedTessellation;

  static ConvexHull convexHull(const std::vector<RealType>& points, 
                               const RealType* low,
                               const RealType& dx) { 
    return ConvexHull(convexHull_3d(points, low, dx), points);
  }
  static IntPoint constructPoint(const RealType* ri,
                                 const RealType* rlow,
                                 const RealType& dx,
                                 const size_t i) {
    return IntPoint(ri[0], ri[1], ri[2], 
                    rlow[0], rlow[1], rlow[2],
                    dx, i);
  }
  static RealPoint constructPoint(const RealType* ri) {
    return RealPoint(ri[0], ri[1], ri[2]);
  }
  static IntPoint faceCentroid(const Tessellation<3, RealType>& mesh,
                               const unsigned iface,
                               const RealType* rlow,
                               const RealType& dx) {
    POLY_ASSERT(iface < mesh.faces.size());
    const unsigned nnodes = mesh.faces[iface].size();
    POLY_ASSERT(nnodes >= 3);
    RealType pface[3] = {0.0, 0.0, 0.0};
    unsigned ni;
    for (typename std::vector<unsigned>::const_iterator itr = mesh.faces[iface].begin();
         itr != mesh.faces[iface].end();
         ++itr) {
      ni = mesh.faces[iface][*itr];
      pface[0] += mesh.nodes[3*ni];
      pface[1] += mesh.nodes[3*ni + 1];
      pface[2] += mesh.nodes[3*ni + 2];
    }
    pface[0] /= nnodes;
    pface[1] /= nnodes;
    pface[2] /= nnodes;
    return constructPoint(pface, rlow, dx, iface);
  }
  static RealType maxLength(const RealType* low, const RealType* high) {
    return std::max(std::max(high[0] - low[0], high[1] - low[1]), high[2] - low[2]);
  }
  static std::vector<RealType> extractCoords(const std::vector<RealType>& allCoords,
                                             const std::vector<unsigned>& indices) {
    std::vector<RealType> result;
    result.reserve(3*indices.size());
    for (std::vector<unsigned>::const_iterator itr = indices.begin();
         itr != indices.end();
         ++itr) {
      const unsigned i = *itr;
      POLY_ASSERT(i < allCoords.size()/3);
      result.push_back(allCoords[3*i]);
      result.push_back(allCoords[3*i + 1]);
      result.push_back(allCoords[3*i + 2]);
    }
    return result;
  }
  static int hullDimension(const ConvexHull hull) {
    POLY_ASSERT(!hull.facets.empty());
    if (hull.facets.size() > 1) {
      int d = 2, i = 0;
      while (d == 2 and i != hull.facets.size()) {
        d = (hull.facets[i].size() == 2) ? 2 : 3;
        ++i;
      }
      return d;
    } else {
      POLY_ASSERT(hull.facets[0].size() == 2);
    return (hull.facets[0][0] == hull.facets[0][1]) ? 0 : 1;
    }
  }
};

}

#endif
