//------------------------------------------------------------------------
// TetgenTessellator
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <limits>
#include <sstream>

#include "polytope.hh" // Pulls in POLY_ASSERT and TetgenTessellator.hh.
#include "Point.hh"
#include "PLC_CSG_3d.hh"
#include "convexHull_2d.hh"
#include "simplifyPLCfacets.hh"
#include "polytope_write_OOGL.hh"
#include "polytope_plc_canned_geometries.hh"

// Pull in tetgen stuff.
#define TETLIBRARY
#include "tetgen.h"

// Returns (positive, 0.0, negative) if pd is (below, coplanar, above) the plane
// (pa, pb, pc), where above is defined such that (pa, pb, pc) is counter-clockwise.
// This puppy is defined in predicates.cc
extern double orient3d(double* pa, double* pb, double* pc, double* pd);

namespace polytope {

using namespace std;
using std::min;
using std::max;
using std::abs;

namespace {

//------------------------------------------------------------------------------
// Emergency dump.
//------------------------------------------------------------------------------
std::string
escapePod(const std::string nameEnd,
          const polytope::PLC<3>& plc,
          const std::vector<double>& points) {
    std::stringstream os;
    os << "test_PLC_" << nameEnd;
    writePLCtoOFF(plc, points, os.str());
    return " : attempted to write to file " + os.str();
}

//------------------------------------------------------------------------------
// Borrow the Point3 type as a tuple to create 3 node facet hashes.
//------------------------------------------------------------------------------
Point3<unsigned>
hashFacet(const unsigned i, const unsigned j, const unsigned k) {
  typedef Point3<unsigned> Tuple3;
  POLY_ASSERT(i != j and i != k and j != k);
  if (i < j and i < k) {
    if (j < k) {
      return Tuple3(i, j, k);
    } else {
      return Tuple3(i, k, j);
    }
  } else if (j < i and j < k) {
    if (i < k) {
      return Tuple3(j, i, k);
    } else {
      return Tuple3(j, k, i);
    }
  } else {
    if (i < j) {
      return Tuple3(k, i, j);
    } else {
      return Tuple3(k, j, i);
    }
  }
}

//------------------------------------------------------------------------------
// Robustly order a convex 3D facet by projecting to 2D, taking the hull, and
// reinserting any collinear edge points along the hull edges.
//------------------------------------------------------------------------------
std::vector<int>
orderFacetNodes3d(const std::vector<int>& facet,
                  const std::vector<double>& points,
                  const double tol) {
  typedef Point2<double> Point2d;
  typedef Point3<double> Point3d;

  std::vector<int> uniqueFacet;
  uniqueFacet.reserve(facet.size());
  std::set<int> used;
  for (std::vector<int>::const_iterator itr = facet.begin();
       itr != facet.end();
       ++itr) {
    if (used.insert(*itr).second) uniqueFacet.push_back(*itr);
  }
  while (uniqueFacet.size() > 1 and uniqueFacet.front() == uniqueFacet.back()) {
    uniqueFacet.pop_back();
  }
  POLY_ASSERT(uniqueFacet.size() >= 3);

  const auto mag2 = [](const Point3d& p) { return geometry::dot<3, double>(&p.x, &p.x); };

  Point3d centroid(0.0, 0.0, 0.0);
  for (unsigned i = 0; i != uniqueFacet.size(); ++i) {
    const unsigned ipoint = uniqueFacet[i];
    centroid.x += points[3*ipoint];
    centroid.y += points[3*ipoint + 1];
    centroid.z += points[3*ipoint + 2];
  }
  centroid /= double(uniqueFacet.size());

  Point3d normal(0.0, 0.0, 0.0);
  for (unsigned i = 1; i + 1 < uniqueFacet.size() and mag2(normal) <= tol*tol; ++i) {
    const unsigned i0 = uniqueFacet[0], i1 = uniqueFacet[i], i2 = uniqueFacet[i + 1];
    const Point3d p0(points[3*i0], points[3*i0 + 1], points[3*i0 + 2]);
    const Point3d p1(points[3*i1], points[3*i1 + 1], points[3*i1 + 2]);
    const Point3d p2(points[3*i2], points[3*i2 + 1], points[3*i2 + 2]);
    const Point3d e1 = p1 - p0;
    const Point3d e2 = p2 - p0;
    geometry::cross<3, double>(&e1.x, &e2.x, &normal.x);
  }
  POLY_ASSERT(mag2(normal) > tol*tol);
  geometry::unitVector<3, double>(&normal.x);

  Point3d axis0(0.0, 0.0, 0.0);
  for (unsigned i = 0; i != uniqueFacet.size() and mag2(axis0) <= tol*tol; ++i) {
    const unsigned ipoint = uniqueFacet[i];
    axis0 = Point3d(points[3*ipoint], points[3*ipoint + 1], points[3*ipoint + 2]) - centroid;
  }
  POLY_ASSERT(mag2(axis0) > tol*tol);
  geometry::unitVector<3, double>(&axis0.x);

  Point3d axis1;
  geometry::cross<3, double>(&normal.x, &axis0.x, &axis1.x);
  POLY_ASSERT(mag2(axis1) > tol*tol);
  geometry::unitVector<3, double>(&axis1.x);

  std::vector<double> coords2d;
  coords2d.reserve(2*uniqueFacet.size());
  double low2d[2] = { std::numeric_limits<double>::max(),  std::numeric_limits<double>::max() };
  for (unsigned i = 0; i != uniqueFacet.size(); ++i) {
    const unsigned ipoint = uniqueFacet[i];
    const Point3d p(points[3*ipoint], points[3*ipoint + 1], points[3*ipoint + 2]);
    const Point3d dp = p - centroid;
    const double u = geometry::dot<3, double>(&dp.x, &axis0.x);
    const double v = geometry::dot<3, double>(&dp.x, &axis1.x);
    coords2d.push_back(u);
    coords2d.push_back(v);
    low2d[0] = std::min(low2d[0], u);
    low2d[1] = std::min(low2d[1], v);
  }

  const PLC<2> hull2d = convexHull_2d(coords2d, low2d, tol);
  POLY_ASSERT(hull2d.facets.size() >= 3);

  std::vector<int> orderedFacet;
  std::vector<unsigned> assigned(uniqueFacet.size(), 0);
  for (unsigned iedge = 0; iedge != hull2d.facets.size(); ++iedge) {
    const unsigned ia = hull2d.facets[iedge][0];
    const unsigned ib = hull2d.facets[iedge][1];
    POLY_ASSERT(ia < uniqueFacet.size());
    POLY_ASSERT(ib < uniqueFacet.size());
    const Point2d a(coords2d[2*ia], coords2d[2*ia + 1]);
    const Point2d b(coords2d[2*ib], coords2d[2*ib + 1]);
    if (orderedFacet.empty()) {
      orderedFacet.push_back(uniqueFacet[ia]);
      assigned[ia] = 1;
    }

    std::vector<std::pair<double, int> > edgePoints;
    const Point2d ab = b - a;
    const double ab2 = std::max(tol*tol, double(geometry::dot<2, double>(&ab.x, &ab.x)));
    for (unsigned j = 0; j != uniqueFacet.size(); ++j) {
      if (j != ia and j != ib and !assigned[j] and
          geometry::between<2, double>(&a.x, &b.x, &coords2d[2*j], tol)) {
        const Point2d p(coords2d[2*j], coords2d[2*j + 1]);
        const Point2d ap = p - a;
        const double t = geometry::dot<2, double>(&ap.x, &ab.x)/ab2;
        edgePoints.push_back(std::make_pair(t, uniqueFacet[j]));
        assigned[j] = 1;
      }
    }
    std::sort(edgePoints.begin(), edgePoints.end());
    for (unsigned j = 0; j != edgePoints.size(); ++j) orderedFacet.push_back(edgePoints[j].second);
    orderedFacet.push_back(uniqueFacet[ib]);
    assigned[ib] = 1;
  }

  POLY_ASSERT(orderedFacet.size() >= 4);
  POLY_ASSERT(orderedFacet.front() == orderedFacet.back());
  orderedFacet.pop_back();
  return orderedFacet;
}

//------------------------------------------------------------------------------
// Given an array of 4 integers and 2 unique values, find the other two.
//------------------------------------------------------------------------------
void
findOtherTetIndices(const int* indices,
                    const int a,
                    const int b,
                    int& c,
                    int& d) {
  POLY_ASSERT(a == indices[0] or a == indices[1] or a == indices[2] or a == indices[3]);
  POLY_ASSERT(b == indices[0] or b == indices[1] or b == indices[2] or b == indices[3]);
  POLY_ASSERT(indices[0] != indices[1] and indices[0] != indices[2] and indices[0] != indices[3] and
                                           indices[1] != indices[2] and indices[1] != indices[3] and
                                                                        indices[2] != indices[3]);
  if (a != indices[0] and b != indices[0]) {
    c = indices[0];
    d = (a != indices[1] and b != indices[1] ? indices[1] :
         a != indices[2] and b != indices[2] ? indices[2] :
         indices[3]);
  } else if (a != indices[1] and b != indices[1]) {
    c = indices[1];
    d = (a != indices[2] and b != indices[2] ? indices[2] :
         indices[3]);
  } else {
    c = indices[2];
    d = indices[3];
  }
}

//------------------------------------------------------------------------------
// Build a ReducedPLC representation of a cell.
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<3, double>
plcOfCell(const internal::QuantTessellation<3, RealType>& qmesh,
          const unsigned icell) {
  typedef typename internal::QuantTessellation<3, RealType>::PointHash PointHash;
  typedef Point3<double> RealPoint;
  ReducedPLC<3, double> result;
  std::map<int, int> old2new;
  for (unsigned i = 0; i != qmesh.cells[icell].size(); ++i) {
    result.facets.push_back(vector<int>());
    if (qmesh.cells[icell][i] < 0) {
      const unsigned iface = ~qmesh.cells[icell][i];
      const unsigned nnodes = qmesh.faces[iface].size();
      for (int j = nnodes - 1; j != -1; --j) {
        const int iedge = qmesh.faces[iface][j];
        const int ip = iedge < 0 ? qmesh.edges[~iedge].first : qmesh.edges[iedge].second;
        if (old2new.find(ip) == old2new.end()) {
          old2new[ip] = result.points.size()/3;
          RealPoint p = qmesh.unhashPosition(qmesh.points[ip]);
          result.points.push_back(p.x);
          result.points.push_back(p.y);
          result.points.push_back(p.z);
        }
        result.facets.back().push_back(old2new[ip]);
      }
      POLY_ASSERT(result.facets.back().size() == nnodes);
    } else {
      const unsigned iface = qmesh.cells[icell][i];
      const unsigned nnodes = qmesh.faces[iface].size();
      for (int j = 0; j != nnodes; ++j) {
        const int iedge = qmesh.faces[iface][j];
        const int ip = iedge < 0 ? qmesh.edges[~iedge].second : qmesh.edges[iedge].first;
        if (old2new.find(ip) == old2new.end()) {
          old2new[ip] = result.points.size()/3;
          RealPoint p = qmesh.unhashPosition(qmesh.points[ip]);
          result.points.push_back(p.x);
          result.points.push_back(p.y);
          result.points.push_back(p.z);
        }
        result.facets.back().push_back(old2new[ip]);
      }
      POLY_ASSERT(result.facets.back().size() == nnodes);
    }
  }
  POLY_ASSERT(result.facets.size() == qmesh.cells[icell].size());
  return result;
}

// //------------------------------------------------------------------------------
// // Compare a point with plane: returns (-1,0,1) for the point 
// // (below,coplanar,above) the plane.
// // For now we require the plane be aligned in the x, y, or z direction: i.e., 
// // normals (+/-1,0,0), (0,+/-1,0), or (0,0,+/-1).
// //------------------------------------------------------------------------------
// template<typename RealType>
// int
// compare(uint64_t point,
//         uint64_t pointPlane_fine,
//         uint64_t pointPlane_coarse,
//         const int64_t xnorm,
//         const int64_t ynorm,
//         const int64_t znorm,
//         const internal::QuantTessellation<3, RealType>& qmesh) {
//   typedef uint64_t PointHash;
//   typedef geometry::Hasher<3, double> HasherType;
//   typedef typename internal::QuantTessellation<3, RealType>::RealPoint RealPoint;
//   const PointHash outerFlag = geometry::Hasher<3, RealType>::outerFlag();
//   POLY_ASSERT(!(pointPlane_fine & outerFlag));
//   POLY_ASSERT(pointPlane_coarse & outerFlag);
//   POLY_ASSERT(std::abs(xnorm) + std::abs(ynorm) + std::abs(znorm) == 1);

//   const uint64_t pointPlane = (point & outerFlag) ? pointPlane_coarse : pointPlane_fine;

//   // Now check that sucker.
//   int testval;
//   if (std::abs(xnorm) == 1) {
//     testval = int64_t(HasherType::qxval(point) - HasherType::qxval(pointPlane))*xnorm;
//   } else if (std::abs(ynorm) == 1) {
//     testval = int64_t(HasherType::qyval(point) - HasherType::qyval(pointPlane))*ynorm;
//   } else {
//     testval = int64_t(HasherType::qzval(point) - HasherType::qzval(pointPlane))*znorm;
//   }
//   return (testval < 0 ? -1 :
//           testval > 0 ?  1 :
//           0);
// }

// //------------------------------------------------------------------------------
// // Same thing for a collection of points.
// //------------------------------------------------------------------------------
// template<typename RealType>
// int
// compare(const std::vector<uint64_t>& points,
//         const uint64_t pointPlane_fine,
//         const uint64_t pointPlane_coarse,
//         const int64_t xnorm,
//         const int64_t ynorm,
//         const int64_t znorm,
//         const internal::QuantTessellation<3, RealType>& qmesh) {
//   typedef uint64_t PointHash;
//   typedef geometry::Hasher<3, double> HasherType;
//   typedef typename internal::QuantTessellation<3, RealType>::RealPoint RealPoint;
//   const PointHash outerFlag = geometry::Hasher<3, RealType>::outerFlag();
//   POLY_ASSERT(!(pointPlane_fine & outerFlag));
//   POLY_ASSERT(pointPlane_coarse & outerFlag);
//   POLY_ASSERT(std::abs(xnorm) + std::abs(ynorm) + std::abs(znorm) == 1);

//   int result = 0, resulti, testval;
//   PointHash pointPlanei;
//   for (unsigned i = 0; i != points.size(); ++i) {

//     // If the test point is quantized in the outer box, we have to translate the 
//     // plane point to the same coarse coordinates.
//     pointPlanei = ((points[i] & outerFlag) ? pointPlane_coarse : pointPlane_fine);

//     // Now check that sucker.
//     if (std::abs(xnorm) == 1) {
//       testval = int64_t(HasherType::qxval(points[i]) - HasherType::qxval(pointPlanei))*xnorm;
//     } else if (std::abs(ynorm) == 1) {
//       testval = int64_t(HasherType::qyval(points[i]) - HasherType::qyval(pointPlanei))*ynorm;
//     } else {
//       testval = int64_t(HasherType::qzval(points[i]) - HasherType::qzval(pointPlanei))*znorm;
//     }
//     resulti = (testval < 0 ? -1 :
//                testval > 0 ?  1 :
//                0);
//     if ((resulti == 0) or
//         (i > 0 and resulti != result)) return 0;
//     result = resulti;
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Quantized plane intersection with a line segment.
// //------------------------------------------------------------------------------
// template<typename RealType>
// uint64_t
// planeIntersection(const uint64_t e1,
//                   const uint64_t e2,
//                   const uint64_t pointPlane,
//                   const int64_t xnorm,
//                   const int64_t ynorm,
//                   const int64_t znorm,
//                   const internal::QuantTessellation<3, RealType>& qmesh) {

//   typedef uint64_t PointHash;
//   typedef geometry::Hasher<3, double> HasherType;
//   typedef typename internal::QuantTessellation<3, RealType>::RealPoint RealPoint;
//   const PointHash outerFlag = geometry::Hasher<3, RealType>::outerFlag();
//   POLY_ASSERT(!(pointPlane & outerFlag));

//   // Export our work to the floating method.
//   // This is definitely not the most efficient thing we could do, but it's expedient
//   // while we get stuff working!
//   // const int64_t e1x = HasherType::qxval(e1), e1y = HasherType::qyval(e1), e1z = HasherType::qzval(e1),
//   //               e2x = HasherType::qxval(e2), e2y = HasherType::qyval(e2), e2z = HasherType::qzval(e2),
//   //               ppx = HasherType::qxval(pointPlane), ppy = HasherType::qyval(pointPlane), ppz = HasherType::qzval(pointPlane);
//   RealPoint p_ray = qmesh.unhashPosition(e1),
//             n_ray = qmesh.unhashPosition(e2) - p_ray,
//             p_plane = qmesh.unhashPosition(pointPlane),
//             n_plane(xnorm, ynorm, znorm),
//             p_intersect;
//   geometry::unitVector<3, RealType>(&n_ray.x);
//   geometry::unitVector<3, RealType>(&n_plane.x);
//   const bool valid = geometry::rayPlaneIntersection(&p_ray.x, &n_ray.x, &p_plane.x, &n_plane.x,
//                                                     qmesh.degeneracy, &p_intersect.x);
//   POLY_ASSERT2(valid, p_ray << " " << n_ray << " " << p_plane << " " << n_plane << " " << p_intersect);

//   // Convert the answer to a quantized coordinate and we're done.
//   return qmesh.hashPosition(p_intersect);
// }

// //------------------------------------------------------------------------------
// // Comparator for sorting points around a face.
// //------------------------------------------------------------------------------
// template<typename RealType>
// struct FacePointComparator {
//   typedef uint64_t PointHash;
//   typedef int64_t CoordHash;
//   typedef typename internal::QuantTessellation<3, RealType>::RealPoint RealPoint;
//   PointHash iorigin;
//   RealPoint origin, normal;
//   const internal::QuantTessellation<3, RealType>& qmesh;
//   FacePointComparator(const PointHash origini,
//                       const CoordHash xnormi,
//                       const CoordHash ynormi,
//                       const CoordHash znormi,
//                       const internal::QuantTessellation<3, RealType>& qmeshi): iorigin(origini),
//                                                                                origin(qmeshi.unhashPosition(origini)),
//                                                                                normal(xnormi, ynormi, znormi),
//                                                                                qmesh(qmeshi) {}
//   bool operator()(const PointHash a, const PointHash b) {  // Looks like operator<
//     if (a == iorigin) {
//       return true;
//     } else if (b == iorigin) {
//       return false;
//     }
//     const RealPoint da = qmesh.unhashPosition(a) - origin,
//                     db = qmesh.unhashPosition(b) - origin;
//     RealPoint da_cross_db;
//     geometry::cross<3, RealType>(&da.x, &db.x, &da_cross_db.x);
//     const RealType test = geometry::dot<3, RealType>(&da_cross_db.x, &normal.x);
//     return test < 0;
//   }
// };

// //------------------------------------------------------------------------------
// // Clip a ReducedPLC with a plane.  The plane is specified in (point, normal)
// // form in the arguments.  
// // For now we require the plane be aligned in the x, y, or z direction: i.e., 
// // normals (+/-1,0,0), (0,+/-1,0), or (0,0,+/-1).
// //------------------------------------------------------------------------------
// template<typename RealType>
// ReducedPLC<3, uint64_t>
// clipReducedPLC(const ReducedPLC<3, uint64_t>& cell,
//                const uint64_t pointPlane_fine,
//                const uint64_t pointPlane_coarse,
//                const int64_t xnorm,
//                const int64_t ynorm,
//                const int64_t znorm,
//                const internal::QuantTessellation<3, RealType>& qmesh) {
//   typedef uint64_t PointHash;
//   typedef geometry::Hasher<3, double> HasherType;
//   const PointHash outerFlag = geometry::Hasher<3, RealType>::outerFlag();
//   POLY_ASSERT(!(pointPlane_fine & outerFlag));
//   POLY_ASSERT(pointPlane_coarse & outerFlag);
//   POLY_ASSERT(std::abs(xnorm) + std::abs(ynorm) + std::abs(znorm) == 1);

//   // Prepare to build the result up.
//   ReducedPLC<3, PointHash> result;
//   // cerr << "--------------------------------------------------------------------------------" << endl
//   //      << "Clipping against " << qmesh.unhashPosition(pointPlane_fine) << " (" << xnorm << " " << ynorm << " " << znorm << ")" << endl;
//   // for (unsigned i = 0; i != cell.facets.size(); ++i) {
//   //   for (unsigned j = 0; j != cell.facets[i].size(); ++j) {
//   //     cerr << "  " << qmesh.unhashPosition(cell.points[cell.facets[i][j]]);
//   //   }
//   //   cerr << endl;
//   // }
//   // cerr << "--------------------------------------------------------------------------------" << endl;

//   // If the entire PLC is on one side of the input plane the answer is simple.
//   const int allcompare = compare(cell.points, pointPlane_fine, pointPlane_coarse, xnorm, ynorm, znorm, qmesh);
//   POLY_ASSERT(allcompare >= -1 and allcompare <= 1);
//   if (allcompare == -1) {
//     return result;
//   } else if (allcompare == 1) {
//     return cell;
//   }

//   // The plane intersects this polyhedron.  Walk each face and intersect it with
//   // the plane.
//   map<PointHash, int> point2id;
//   vector<PointHash> newfacet;  // We may create at most one new facet by clipping with a plane.
//   for (unsigned iface = 0; iface != cell.facets.size(); ++iface) {
//     const unsigned nnodes = cell.facets[iface].size();
//     POLY_ASSERT(nnodes >= 3);
//     vector<int> clippedfacet;
//     for (unsigned i = 0; i != nnodes; ++i) {
//       const unsigned e1 = cell.facets[iface][i],
//                      e2 = cell.facets[iface][(i + 1) % nnodes];
//       const int compare1 = compare(cell.points[e1], pointPlane_fine, pointPlane_coarse, xnorm, ynorm, znorm, qmesh),
//                 compare2 = compare(cell.points[e2], pointPlane_fine, pointPlane_coarse, xnorm, ynorm, znorm, qmesh);

//       // Does the first point of this edge make the cut?
//       if (compare1 != -1) {
//         const unsigned oldsize = point2id.size();
//         const unsigned k = internal::addKeyToMap(cell.points[e1], point2id);
//         if (k == oldsize) result.points.push_back(cell.points[e1]);
//         if ((clippedfacet.size() == 0) or
//             ((i < nnodes - 1) and (k != clippedfacet.back())) or
//             ((i == nnodes - 1) and (k != clippedfacet.back()) and (k != clippedfacet[0]))) clippedfacet.push_back(k);
//       }

//       // Does this edge stradle the plane?
//       if (compare1*compare2 == -1) {
//         const PointHash newpoint = planeIntersection(cell.points[e1],
//                                                      cell.points[e2],
//                                                      pointPlane_fine,
//                                                      xnorm,
//                                                      ynorm,
//                                                      znorm,
//                                                      qmesh);
//         const unsigned oldsize = point2id.size();
//         const unsigned k = internal::addKeyToMap(newpoint, point2id);
//         if (k == oldsize) result.points.push_back(newpoint);
//         if ((clippedfacet.size() == 0) or
//             ((i < nnodes - 1) and (k != clippedfacet.back())) or
//             ((i == nnodes - 1) and (k != clippedfacet.back()) and (k != clippedfacet[0]))) clippedfacet.push_back(k);
//         newfacet.push_back(newpoint);
//       }
//     }

//     // Make sure there are no repeats in the clipped facet.
//     if (clippedfacet.size() >= 3) {
//       vector<int> clippedfacet2;
//       for (int i = 0; i != clippedfacet.size(); ++i) {
//         int j = (i + 1) % clippedfacet.size();
//         if (clippedfacet[i] != clippedfacet[j]) clippedfacet2.push_back(clippedfacet[i]);
//       }
//       if (clippedfacet2.size() >= 3) result.facets.push_back(clippedfacet2);
//     }
//   }

//   // Did we create a new facet?
//   if (newfacet.size() >= 3) {
//     std::sort(newfacet.begin(), newfacet.end(), FacePointComparator<RealType>(newfacet.front(), xnorm, ynorm, znorm, qmesh));
//     newfacet.erase(std::unique(newfacet.begin(), newfacet.end()), newfacet.end());
//     if (newfacet.size() >= 3) {
//       // cerr << "New facet points : ";
//       // for (unsigned i = 0; i != newfacet.size(); ++i) cerr << " " << point2id[newfacet[i]];
//       // cerr << endl;
//       // cerr << "New facet points ordering : ";
//       // FacePointComparator<RealType> comp(newfacet.front(), xnorm, ynorm, znorm, qmesh);
//       // for (unsigned i = 0; i != newfacet.size(); ++i) cerr << " "
//       //                                                      << comp(newfacet[i], newfacet[(i + 1) % newfacet.size()]);
//       // cerr << endl;
//       vector<int> newfacetIDs;
//       for (unsigned i = 0; i != newfacet.size(); ++i) newfacetIDs.push_back(point2id[newfacet[i]]);
//       result.facets.push_back(newfacetIDs);
//     }
//   }

//   // That's it.
//   POLY_ASSERT(result.facets.size() >= 4);
//   return result;
// }

} // end anonymous namespace

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
TetgenTessellator::
TetgenTessellator():
  Tessellator<3, double>() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
TetgenTessellator::
~TetgenTessellator() {
}

//------------------------------------------------------------------------------
// Quantized tessellation.
//------------------------------------------------------------------------------
void
TetgenTessellator::
tessellateQuantized(QuantizedTessellation& /*qmesh*/) const {
  POLY_ASSERT2(false, "TetgenTessellator has not been migrated to the QuantizedTessellation API yet. Use the direct 3D tessellate overloads for now.");
}

//------------------------------------------------------------------------------
// Unbounded tessellation.
//------------------------------------------------------------------------------
void
TetgenTessellator::
tessellate(const vector<double>& points,
           Tessellation<3, double>& mesh) const {

  // First generate our internal quantized tessellation representation.
  internal::QuantTessellation<3, double> qmesh;
  vector<double> nonGeneratingPoints;
  this->computeUnboundedQuantizedTessellation(points, nonGeneratingPoints, qmesh);

  // Convert to the output tessellation and we're done.
  qmesh.tessellation(mesh);
}

//------------------------------------------------------------------------------
// Tessellate within a bounding box.
//------------------------------------------------------------------------------
void
TetgenTessellator::
tessellate(const vector<double>& points,
           double* low,
           double* high,
           Tessellation<3, double>& mesh) const {

  // Create a PLC of the bounding points, and use the ReducedPLC method to do the 
  // tessellation.
  ReducedPLC<3, double> box = plc_box<3, double>(low, high);
  this->tessellate(points, box, mesh);

  // Make a final pass and ensure our bounds are met.
  for (unsigned i = 0; i != mesh.nodes.size()/3; ++i) {
    mesh.nodes[3*i  ] = std::max(low[0], std::min(high[0], mesh.nodes[3*i  ]));
    mesh.nodes[3*i+1] = std::max(low[1], std::min(high[1], mesh.nodes[3*i+1]));
    mesh.nodes[3*i+2] = std::max(low[2], std::min(high[2], mesh.nodes[3*i+2]));
  }
}

//------------------------------------------------------------------------------
// Tessellate within a PLC.
//------------------------------------------------------------------------------
void
TetgenTessellator::
tessellate(const vector<double>& points,
           const vector<double>& PLCpoints,
           const PLC<3>& geometry,
           Tessellation<3, double>& mesh) const {

  // We export to the ReducedPLC method.
  ReducedPLC<3, double> boundary;
  boundary.facets = geometry.facets;
  boundary.holes = geometry.holes;
  boundary.points = PLCpoints;
  this->tessellate(points, boundary, mesh);
}

//------------------------------------------------------------------------------
// Tessellate within a ReducedPLC.
//------------------------------------------------------------------------------
void
TetgenTessellator::
tessellate(const vector<double>& points,
           const ReducedPLC<3, double>& geometry,
           Tessellation<3, double>& mesh) const {

  typedef geometry::Hasher<3, double> HasherType;
  typedef internal::QuantTessellation<3, double>::PointHash PointHash;
  typedef internal::QuantTessellation<3, double>::EdgeHash EdgeHash;
  typedef internal::QuantTessellation<3, double>::IntPoint IntPoint;
  typedef internal::QuantTessellation<3, double>::RealPoint RealPoint;
  typedef geometry::Hasher<3, double> HasherType;

  // const PointHash outerFlag = HasherType::outerFlag();

  escapePod("boundary", geometry, geometry.points);

  // Pre-conditions.
  POLY_ASSERT(not points.empty());
  POLY_ASSERT(points.size() % 3 == 0);
  const unsigned numGenerators = points.size()/3;
  // for (unsigned i = 0; i != numGenerators; ++i) {
  //   for (unsigned j = 0; j != 3; ++j) {
  //     POLY_ASSERT(low[j] <= points[3*i+j] and points[3*i+j] <= high[j]);
  //   }
  // }

  // Create the unbounded QuantTessellation.
  internal::QuantTessellation<3, double> qmesh0;
  this->computeUnboundedQuantizedTessellation(points, geometry.points, qmesh0);

  // Create a new QuantTessellation.  This one will only use the single level of
  // quantization since we know the PLC is within this inner region.
  internal::QuantTessellation<3, double> qmesh1;
  qmesh1.generators = qmesh0.generators;
  qmesh1.low_labframe = qmesh0.low_labframe;
  qmesh1.high_labframe = qmesh0.high_labframe;
  qmesh1.low_inner = qmesh0.low_inner;
  qmesh1.high_inner = qmesh0.high_inner;
  qmesh1.low_outer = qmesh0.low_inner;
  qmesh1.high_outer = qmesh0.high_inner;
  qmesh1.degeneracy = 1.0e-5;

  // Walk each of the cells in the unbounded tessellation.
  for (unsigned icell = 0; icell != numGenerators; ++icell) {
    // Build a PLC to represent just this cell.
    ReducedPLC<3, double> cell = plcOfCell(qmesh0, icell);

    // Intersect with the boundary.
    cell = CSG::csg_intersect(geometry, cell);
    cell = simplifyPLCfacets(cell, cell.points, &qmesh1.low_inner.x, &qmesh1.high_inner.x, 1.0e-5);   // We have to use a *much* coarser degeneracy here due to CSG accuracy...  :(
    POLY_ASSERT(cell.facets.size() >= 4);

    // Add this cell and its elements to the new tessellation.
    vector<int> nodeIDs, edgeIDs, faceIDs;
    qmesh1.cells.push_back(vector<int>());
    for (unsigned i = 0; i != cell.points.size(); ++i) {
      nodeIDs.push_back(qmesh1.addNewNode(HasherType::hashPosition(&cell.points[3*i],
                                                                   const_cast<double*>(&qmesh1.low_inner.x), const_cast<double*>(&qmesh1.high_inner.x),
                                                                   const_cast<double*>(&qmesh1.low_outer.x), const_cast<double*>(&qmesh1.high_outer.x),
                                                                   qmesh1.degeneracy)));
    }
    for (unsigned iface = 0; iface != cell.facets.size(); ++iface) {
      vector<int> facetNodes = orderFacetNodes3d(cell.facets[iface], cell.points, qmesh1.degeneracy);

      // Remove vertices that lie on straight runs of the polygon boundary.
      bool removedCollinear = true;
      while (removedCollinear and facetNodes.size() > 3) {
        removedCollinear = false;
        vector<int> cleaned;
        cleaned.reserve(facetNodes.size());
        for (unsigned i = 0; i != facetNodes.size(); ++i) {
          const unsigned iprev = facetNodes[(i + facetNodes.size() - 1) % facetNodes.size()];
          const unsigned icurr = facetNodes[i];
          const unsigned inext = facetNodes[(i + 1) % facetNodes.size()];
          const RealPoint pprev(cell.points[3*iprev], cell.points[3*iprev + 1], cell.points[3*iprev + 2]);
          const RealPoint pcurr(cell.points[3*icurr], cell.points[3*icurr + 1], cell.points[3*icurr + 2]);
          const RealPoint pnext(cell.points[3*inext], cell.points[3*inext + 1], cell.points[3*inext + 2]);
          const RealPoint eprev = pcurr - pprev;
          const RealPoint enext = pnext - pcurr;
          RealPoint turn;
          geometry::cross<3, double>(&eprev.x, &enext.x, &turn.x);
          if (geometry::dot<3, double>(&turn.x, &turn.x) > 1.0e-20) {
            cleaned.push_back(icurr);
          } else {
            removedCollinear = true;
          }
        }
        facetNodes.swap(cleaned);
      }
      POLY_ASSERT(facetNodes.size() >= 3);

      // Orient the local facet cycle so this cell sees a positive signed
      // volume. Shared faces reused by other cells can then derive the proper
      // opposite sign from the stored face orientation.
      {
        const RealPoint p0(cell.points[3*facetNodes[0]], cell.points[3*facetNodes[0] + 1], cell.points[3*facetNodes[0] + 2]);
        const RealPoint p1(cell.points[3*facetNodes[1]], cell.points[3*facetNodes[1] + 1], cell.points[3*facetNodes[1] + 2]);
        const RealPoint p2(cell.points[3*facetNodes[2]], cell.points[3*facetNodes[2] + 1], cell.points[3*facetNodes[2] + 2]);
        const RealType localVol = geometry::tetrahedralVolume6(&qmesh1.generators[3*icell], &p2.x, &p1.x, &p0.x);
        POLY_ASSERT(localVol != 0.0);
        if (localVol < 0.0) std::reverse(facetNodes.begin(), facetNodes.end());
      }

      vector<int> face;
      face.reserve(facetNodes.size());
      for (unsigned i = 0; i != facetNodes.size(); ++i) {
        const unsigned j = (i + 1) % facetNodes.size();
        const EdgeHash ehash = internal::hashEdge(nodeIDs[facetNodes[i]],
                                                  nodeIDs[facetNodes[j]]);
        int iedge = qmesh1.addNewEdge(ehash);
        if (ehash.first == nodeIDs[facetNodes[j]]) iedge = ~iedge;
        face.push_back(iedge);
      }
      POLY_ASSERT(face.size() == facetNodes.size());
      const unsigned meshFace = qmesh1.addNewFace(face);
      POLY_ASSERT(meshFace < qmesh1.faces.size());

      // Choose the face sign from the stored face orientation, not from whether
      // the face was newly inserted. Shared bounded faces are reused by multiple
      // cells, and the second cell needs the opposite sign if the stored face
      // normal points away from its generator.
      POLY_ASSERT(qmesh1.faces[meshFace].size() >= 3);
      std::vector<unsigned> faceNodes;
      faceNodes.reserve(qmesh1.faces[meshFace].size());
      for (unsigned j = 0; j != qmesh1.faces[meshFace].size(); ++j) {
        const int iedge = qmesh1.faces[meshFace][j];
        faceNodes.push_back(iedge >= 0 ? qmesh1.edges[iedge].first : qmesh1.edges[~iedge].second);
      }
      POLY_ASSERT(faceNodes.size() >= 3);
      const RealPoint p0 = qmesh1.nodePosition(faceNodes[0]);
      const RealPoint p1 = qmesh1.nodePosition(faceNodes[1]);
      const RealPoint p2 = qmesh1.nodePosition(faceNodes[2]);
      const RealType vol = geometry::tetrahedralVolume6(&qmesh1.generators[3*icell], &p2.x, &p1.x, &p0.x);
      POLY_ASSERT(vol != 0.0);
      qmesh1.cells.back().push_back(vol > 0.0 ? meshFace : ~meshFace);
    }
    POLY_ASSERT(qmesh1.cells.back().size() == cell.facets.size());
  }

  // The QuantTessellation should be complete now.
  qmesh1.assertValid();

  // Convert to the output tessellation and we're done.
  qmesh1.tessellation(mesh);
}

// Alternate version of ReducedPLC tessellation
// void
// TetgenTessellator::
// tessellate(const vector<double>& points,
//            const ReducedPLC<3, double>& geometry,
//            Tessellation<3, double>& mesh) const {

//   typedef geometry::Hasher<3, double> HasherType;
//   typedef internal::QuantTessellation<3, double>::PointHash PointHash;
//   typedef internal::QuantTessellation<3, double>::EdgeHash EdgeHash;
//   typedef internal::QuantTessellation<3, double>::IntPoint IntPoint;
//   typedef internal::QuantTessellation<3, double>::RealPoint RealPoint;
//   typedef geometry::Hasher<3, double> HasherType;

//   // const PointHash outerFlag = HasherType::outerFlag();

//   escapePod("boundary", geometry, geometry.points);

//   // Pre-conditions.
//   POLY_ASSERT(not points.empty());
//   POLY_ASSERT(points.size() % 3 == 0);
//   const unsigned numGenerators = points.size()/3;
//   // for (unsigned i = 0; i != numGenerators; ++i) {
//   //   for (unsigned j = 0; j != 3; ++j) {
//   //     POLY_ASSERT(low[j] <= points[3*i+j] and points[3*i+j] <= high[j]);
//   //   }
//   // }

//   // Create the unbounded QuantTessellation.
//   internal::QuantTessellation<3, double> qmesh0;
//   this->computeUnboundedQuantizedTessellation(points, geometry.points, qmesh0);

//   // Create a new QuantTessellation.  This one will only use the single level of
//   // quantization since we know the PLC is within this inner region.
//   internal::QuantTessellation<3, double> qmesh1;
//   qmesh1.generators = qmesh0.generators;
//   qmesh1.low_labframe = qmesh0.low_labframe;
//   qmesh1.high_labframe = qmesh0.high_labframe;
//   qmesh1.low_inner = qmesh0.low_inner;
//   qmesh1.high_inner = qmesh0.high_inner;
//   qmesh1.low_outer = qmesh0.low_inner;
//   qmesh1.high_outer = qmesh0.high_inner;
//   qmesh1.degeneracy = 1.0e-5;

//   // Walk each of the cells in the unbounded tessellation.
//   for (unsigned icell = 0; icell != numGenerators; ++icell) {
//     // Build a PLC to represent just this cell.
//     ReducedPLC<3, double> cell = plcOfCell(qmesh0, icell);

//     // Intersect with the boundary.
//     cell = CSG::csg_intersect(geometry, cell);
//     cell = simplifyPLCfacets(cell, cell.points, &qmesh1.low_inner.x, &qmesh1.high_inner.x, 1.0e-5);   // We have to use a *much* coarser degeneracy here due to CSG accuracy...  :(
//     POLY_ASSERT(cell.facets.size() >= 4);

//     // Add this cell and its elements to the new tessellation.
//     vector<int> nodeIDs, edgeIDs, faceIDs;
//     qmesh1.cells.push_back(vector<int>());
//     for (unsigned i = 0; i != cell.points.size(); ++i) {
//       nodeIDs.push_back(qmesh1.addNewNode(HasherType::hashPosition(&cell.points[3*i],
//                                                                    const_cast<double*>(&qmesh1.low_inner.x), const_cast<double*>(&qmesh1.high_inner.x),
//                                                                    const_cast<double*>(&qmesh1.low_outer.x), const_cast<double*>(&qmesh1.high_outer.x),
//                                                                    qmesh1.degeneracy)));
//     }
//     for (unsigned iface = 0; iface != cell.facets.size(); ++iface) {
//       vector<int> facetNodes = orderFacetNodes3d(cell.facets[iface], cell.points, qmesh1.degeneracy);

//       // Remove vertices that lie on straight runs of the polygon boundary.
//       bool removedCollinear = true;
//       while (removedCollinear and facetNodes.size() > 3) {
//         removedCollinear = false;
//         vector<int> cleaned;
//         cleaned.reserve(facetNodes.size());
//         for (unsigned i = 0; i != facetNodes.size(); ++i) {
//           const unsigned iprev = facetNodes[(i + facetNodes.size() - 1) % facetNodes.size()];
//           const unsigned icurr = facetNodes[i];
//           const unsigned inext = facetNodes[(i + 1) % facetNodes.size()];
//           const RealPoint pprev(cell.points[3*iprev], cell.points[3*iprev + 1], cell.points[3*iprev + 2]);
//           const RealPoint pcurr(cell.points[3*icurr], cell.points[3*icurr + 1], cell.points[3*icurr + 2]);
//           const RealPoint pnext(cell.points[3*inext], cell.points[3*inext + 1], cell.points[3*inext + 2]);
//           const RealPoint eprev = pcurr - pprev;
//           const RealPoint enext = pnext - pcurr;
//           RealPoint turn;
//           geometry::cross<3, double>(&eprev.x, &enext.x, &turn.x);
//           if (geometry::dot<3, double>(&turn.x, &turn.x) > 1.0e-20) {
//             cleaned.push_back(icurr);
//           } else {
//             removedCollinear = true;
//           }
//         }
//         facetNodes.swap(cleaned);
//       }
//       POLY_ASSERT(facetNodes.size() >= 3);

//       // Orient the local facet cycle so this cell sees a positive signed
//       // volume. Shared faces reused by other cells can then derive the proper
//       // opposite sign from the stored face orientation.
//       {
//         const RealPoint p0(cell.points[3*facetNodes[0]], cell.points[3*facetNodes[0] + 1], cell.points[3*facetNodes[0] + 2]);
//         const RealPoint p1(cell.points[3*facetNodes[1]], cell.points[3*facetNodes[1] + 1], cell.points[3*facetNodes[1] + 2]);
//         const RealPoint p2(cell.points[3*facetNodes[2]], cell.points[3*facetNodes[2] + 1], cell.points[3*facetNodes[2] + 2]);
//         const RealType localVol = geometry::tetrahedralVolume6(&qmesh1.generators[3*icell], &p2.x, &p1.x, &p0.x);
//         POLY_ASSERT(localVol != 0.0);
//         if (localVol < 0.0) std::reverse(facetNodes.begin(), facetNodes.end());
//       }

//       vector<int> face;
//       face.reserve(facetNodes.size());
//       for (unsigned i = 0; i != facetNodes.size(); ++i) {
//         const unsigned j = (i + 1) % facetNodes.size();
//         const EdgeHash ehash = internal::hashEdge(nodeIDs[facetNodes[i]],
//                                                   nodeIDs[facetNodes[j]]);
//         int iedge = qmesh1.addNewEdge(ehash);
//         if (ehash.first == nodeIDs[facetNodes[j]]) iedge = ~iedge;
//         face.push_back(iedge);
//       }
//       POLY_ASSERT(face.size() == facetNodes.size());
//       const unsigned meshFace = qmesh1.addNewFace(face);
//       POLY_ASSERT(meshFace < qmesh1.faces.size());

//       // Choose the face sign from the stored face orientation, not from whether
//       // the face was newly inserted. Shared bounded faces are reused by multiple
//       // cells, and the second cell needs the opposite sign if the stored face
//       // normal points away from its generator.
//       POLY_ASSERT(qmesh1.faces[meshFace].size() >= 3);
//       std::vector<unsigned> faceNodes;
//       faceNodes.reserve(qmesh1.faces[meshFace].size());
//       for (unsigned j = 0; j != qmesh1.faces[meshFace].size(); ++j) {
//         const int iedge = qmesh1.faces[meshFace][j];
//         faceNodes.push_back(iedge >= 0 ? qmesh1.edges[iedge].first : qmesh1.edges[~iedge].second);
//       }
//       POLY_ASSERT(faceNodes.size() >= 3);
//       const RealPoint p0 = qmesh1.nodePosition(faceNodes[0]);
//       const RealPoint p1 = qmesh1.nodePosition(faceNodes[1]);
//       const RealPoint p2 = qmesh1.nodePosition(faceNodes[2]);
//       const RealType vol = geometry::tetrahedralVolume6(&qmesh1.generators[3*icell], &p2.x, &p1.x, &p0.x);
//       POLY_ASSERT(vol != 0.0);
//       qmesh1.cells.back().push_back(vol > 0.0 ? meshFace : ~meshFace);
//     }
//     POLY_ASSERT(qmesh1.cells.back().size() == cell.facets.size());
//   }

//   // The QuantTessellation should be complete now.
//   qmesh1.assertValid();

//   // Convert to the output tessellation and we're done.
//   qmesh1.tessellation(mesh);
// }

//------------------------------------------------------------------------------
// Internal method that returns an intermediated quantized representation
// of the unbounded tessellation.
//------------------------------------------------------------------------------
void
TetgenTessellator::
computeUnboundedQuantizedTessellation(const vector<double>& points,
                                      const vector<double>& nonGeneratingPoints,
                                      internal::QuantTessellation<3, double>& qmesh) const {

  // Pre-conditions.
  POLY_ASSERT(not points.empty());
  POLY_ASSERT(points.size() % 3 == 0);
  POLY_ASSERT(nonGeneratingPoints.size() % 3 == 0);

  typedef internal::QuantTessellation<3, double>::PointHash PointHash;
  typedef internal::QuantTessellation<3, double>::EdgeHash EdgeHash;
  typedef Point3<double> RealPoint;
  typedef Point3<unsigned> TetFacetHash;  // kind of nefarious!

  qmesh.degeneracy = mDegeneracy;

  // Compute the normalized generators.
  const unsigned numGenerators = points.size() / 3;
  //qmesh.generators = this->computeNormalizedPoints(points, nonGeneratingPoints, true, &qmesh.low_labframe.x, &qmesh.high_labframe.x);
  unsigned i, j, k;

  // Build the input to tetgen.
  tetgenio in;
  in.firstnumber = 0;
  in.mesh_dim = 3;
  in.pointlist = new double[qmesh.generators.size()];
  copy(&qmesh.generators.front(), &qmesh.generators.front() + qmesh.generators.size(), in.pointlist);
  in.pointattributelist = 0;
  in.pointmtrlist = 0;
  in.pointmarkerlist = 0;
  in.numberofpoints = qmesh.generators.size() / 3;
  in.numberofpointattributes = 0;
  in.numberofpointmtrs = 0;

  // Do the tetrahedralization.
  tetgenio out;
  tetrahedralize((char*)"Q", &in, &out);
  // tetrahedralize((char*)"V", &in, &out);

  // Make sure we got something.
  if (out.numberoftetrahedra == 0)
    error("TetgenTessellator: Delauney tetrahedralization produced 0 tetrahedra!");
//  if (out.numberofpoints != numGenerators) {
//    char err[1024];
//    snprintf(err, 1024, "TetgenTessellator: Delauney tetrahedralization produced %d tetrahedra\n(%d generating points given)", 
//             out.numberofpoints, (int)numGenerators);
//    error(err);
//  }

  // Compute the circumcenters of the tetrahedra, and the set of tets associated
  // with each generator.
  qmesh.low_inner = RealPoint(0, 0, 0);
  qmesh.high_inner = RealPoint(1, 1, 1);
  qmesh.low_outer = RealPoint(numeric_limits<double>::max(),
                              numeric_limits<double>::max(),
                              numeric_limits<double>::max());
  qmesh.high_outer = RealPoint(-numeric_limits<double>::max(),
                               -numeric_limits<double>::max(),
                               -numeric_limits<double>::max());
  vector<RealPoint> circumcenters(out.numberoftetrahedra);
  int a, b, c, d;
  EdgeHash ab, ac, ad, bc, bd, cd;
  TetFacetHash abc, abd, bcd, acd;
  map<TetFacetHash, vector<unsigned> > facet2tets;      // Tets which share a facet.
  map<EdgeHash, set<unsigned> > edge2tets;              // Tets which share a given edge.
  for (i = 0; i != out.numberoftetrahedra; ++i) {
    a = out.tetrahedronlist[4*i];
    b = out.tetrahedronlist[4*i+1];
    c = out.tetrahedronlist[4*i+2];
    d = out.tetrahedronlist[4*i+3];
    POLY_ASSERT(a < numGenerators);
    POLY_ASSERT(b < numGenerators);
    POLY_ASSERT(c < numGenerators);
    POLY_ASSERT(d < numGenerators);
    geometry::computeCircumcenter3d(&out.pointlist[3*a],
                                    &out.pointlist[3*b],
                                    &out.pointlist[3*c],
                                    &out.pointlist[3*d],
                                    &circumcenters[i].x);
    ab = internal::hashEdge(a, b);
    ac = internal::hashEdge(a, c);
    ad = internal::hashEdge(a, d);
    bc = internal::hashEdge(b, c);
    bd = internal::hashEdge(b, d);
    cd = internal::hashEdge(c, d);
    abc = hashFacet(a, b, c);
    abd = hashFacet(a, b, d);
    bcd = hashFacet(b, c, d);
    acd = hashFacet(a, c, d);
    facet2tets[abc].push_back(i);
    facet2tets[abd].push_back(i);
    facet2tets[bcd].push_back(i);
    facet2tets[acd].push_back(i);
    edge2tets[ab].insert(i);
    edge2tets[ac].insert(i);
    edge2tets[ad].insert(i);
    edge2tets[bc].insert(i);
    edge2tets[bd].insert(i);
    edge2tets[cd].insert(i);
    qmesh.low_outer.x = min(qmesh.low_outer.x, circumcenters[i].x);
    qmesh.low_outer.y = min(qmesh.low_outer.y, circumcenters[i].y);
    qmesh.low_outer.z = min(qmesh.low_outer.z, circumcenters[i].z);
    qmesh.high_outer.x = max(qmesh.high_outer.x, circumcenters[i].x);
    qmesh.high_outer.y = max(qmesh.high_outer.y, circumcenters[i].y);
    qmesh.high_outer.z = max(qmesh.high_outer.z, circumcenters[i].z);
  }
  POLY_BEGIN_CONTRACT_SCOPE;
  {
    for (map<TetFacetHash, vector<unsigned> >::const_iterator itr = facet2tets.begin();
         itr != facet2tets.end();
         ++itr) POLY_ASSERT(itr->second.size() == 1 or itr->second.size() == 2);
    for (map<EdgeHash, set<unsigned> >::const_iterator itr = edge2tets.begin();
         itr != edge2tets.end();
         ++itr) POLY_ASSERT(itr->second.size() >= 1);
    POLY_ASSERT(qmesh.low_outer.x < qmesh.high_outer.x and
                qmesh.low_outer.y < qmesh.high_outer.y and
                qmesh.low_outer.z < qmesh.high_outer.z);
  }
  POLY_END_CONTRACT_SCOPE;

  // Expand the outer bounding box, and choose our infSphere radius.
  qmesh.low_outer.x = min(qmesh.low_outer.x, qmesh.low_inner.x);
  qmesh.low_outer.y = min(qmesh.low_outer.y, qmesh.low_inner.y);
  qmesh.low_outer.z = min(qmesh.low_outer.z, qmesh.low_inner.z);
  qmesh.high_outer.x = max(qmesh.high_outer.x, qmesh.high_inner.x);
  qmesh.high_outer.y = max(qmesh.high_outer.y, qmesh.high_inner.y);
  qmesh.high_outer.z = max(qmesh.high_outer.z, qmesh.high_inner.z);
  double rinf = 4.0*max(    qmesh.high_outer.x - qmesh.low_outer.x,
                          max(qmesh.high_outer.y - qmesh.low_outer.y,
                              qmesh.high_outer.z - qmesh.low_outer.z));
  const RealPoint centroid_outer = (qmesh.low_outer + qmesh.high_outer)/2;
  qmesh.low_outer.x = centroid_outer.x - 1.05*rinf;
  qmesh.low_outer.y = centroid_outer.y - 1.05*rinf;
  qmesh.low_outer.z = centroid_outer.z - 1.05*rinf;
  qmesh.high_outer.x = centroid_outer.x + 1.05*rinf;
  qmesh.high_outer.y = centroid_outer.y + 1.05*rinf;
  qmesh.high_outer.z = centroid_outer.z + 1.05*rinf;

  // Create the quantized circumcenters, and the map from the (possibly) degenerate
  // circumcenters to their unique IDs.
  map<int, unsigned> tet2id;
  for (i = 0; i != out.numberoftetrahedra; ++i) {
    tet2id[i] = qmesh.addNewNode(circumcenters[i]);
  }

  // Any surface facets create new "infinite" or "unbounded" rays, which originate at
  // the tet circumcenter and pass through the circumcenter of the triangular facet.
  // Look for any surface facets we need to project unbounded rays through.
  bool test;
  RealPoint fhat, tetcent, test_point, a_b, a_c, pinf;
  map<TetFacetHash, unsigned> facet2id;
  qmesh.infNodes = vector<unsigned>();
  for (map<TetFacetHash, vector<unsigned> >::const_iterator facetItr = facet2tets.begin();
       facetItr != facet2tets.end();
       ++facetItr) {
    const TetFacetHash& facet = facetItr->first;
    const vector<unsigned>& tets = facetItr->second;
    if (tets.size() == 1) {
      i = tets[0];
      POLY_ASSERT(i < out.numberoftetrahedra);
      a = out.tetrahedronlist[4*i];
      b = out.tetrahedronlist[4*i+1];
      c = out.tetrahedronlist[4*i+2];
      d = out.tetrahedronlist[4*i+3];
      POLY_ASSERT(a < numGenerators);
      POLY_ASSERT(b < numGenerators);
      POLY_ASSERT(c < numGenerators);
      POLY_ASSERT(d < numGenerators);
      geometry::computeTetCentroid(&out.pointlist[3*a],
                                   &out.pointlist[3*b],
                                   &out.pointlist[3*c],
                                   &out.pointlist[3*d],
                                   &tetcent.x);

      // We need the ray unit vector.
      test = geometry::computeTriangleCircumcenter3d(&out.pointlist[3*facet.x],
                                                     &out.pointlist[3*facet.y],
                                                     &out.pointlist[3*facet.z],
                                                     &fhat.x);
      POLY_ASSERT(test);
      fhat -= circumcenters[i];

      // Check for the special case of the tet circumcenter coplanar with the facet.
      if (abs(geometry::dot<3, RealType>(&fhat.x, &fhat.x)) < mDegeneracy) {
        // Yep, it's in the plane.  Just project the ray out orthogonally to the facet.
        a_b.x = out.pointlist[3*facet.y]   - out.pointlist[3*facet.x];
        a_b.y = out.pointlist[3*facet.y+1] - out.pointlist[3*facet.x+1];
        a_b.z = out.pointlist[3*facet.y+2] - out.pointlist[3*facet.x+2];
        a_c.x = out.pointlist[3*facet.z]   - out.pointlist[3*facet.x];
        a_c.y = out.pointlist[3*facet.z+1] - out.pointlist[3*facet.x+1];
        a_c.z = out.pointlist[3*facet.z+2] - out.pointlist[3*facet.x+2];
        geometry::cross<3, RealType>(&a_b.x, &a_c.x, &fhat.x);
      }
      geometry::unitVector<3, RealType>(&fhat.x);

      // The ray unit vector should point in the opposite direction from the facet as the tet centroid.
      POLY_ASSERT(abs(orient3d(&out.pointlist[3*facet.x], &out.pointlist[3*facet.y], &out.pointlist[3*facet.z], &tetcent.x)) > mDegeneracy);
      copy(&out.pointlist[3*facet.x], &out.pointlist[3*facet.x] + 3, &test_point.x);
      test_point += fhat*rinf;
      if (orient3d(&out.pointlist[3*facet.x], &out.pointlist[3*facet.y], &out.pointlist[3*facet.z], &tetcent.x)*
          orient3d(&out.pointlist[3*facet.x], &out.pointlist[3*facet.y], &out.pointlist[3*facet.z], &test_point.x) > 0.0) fhat *= -1.0;

      // Now we can compute the point where this ray intersects the surrounding "inf" sphere.
      test = geometry::raySphereIntersection(&circumcenters[i].x,
                                             &fhat.x,
                                             &centroid_outer.x,
                                             rinf,
                                             1.0e-10,
                                             &pinf.x);
      POLY_ASSERT(test);
      
      // Add this infPoint to the quantized tessellation.
      k = qmesh.point2id.size();
      j = qmesh.addNewNode(pinf);
      POLY_ASSERT(facet2id.find(facet) == facet2id.end());
      facet2id[facet] = j;
      if (k != qmesh.point2id.size()) qmesh.infNodes.push_back(j);
    }
  }

  // Build the edges and faces corresponding to each tet edge.  Recall here that a tet edge is 
  // actualy the line connecting two generators, so not the edge of the mesh we want.
  int iedge, iface;
  RealPoint ghat, e0, e1, e2, f1, f2;
  RealType vol;
  map<EdgeHash, int> faceMap;
  qmesh.faces.reserve(edge2tets.size());
  qmesh.cells = vector<vector<int> >(numGenerators);
  TetFacetHash lastFacet;
  unsigned ii, jj;
  vector<vector<EdgeHash> > cellInfEdges(numGenerators);
  for (map<EdgeHash, set<unsigned> >::const_iterator edgeItr = edge2tets.begin();
       edgeItr != edge2tets.end();
       ++edgeItr) {
    const EdgeHash& ehash = edgeItr->first;
    const set<unsigned>& tets = edgeItr->second;
    a = ehash.first;
    b = ehash.second;
    POLY_ASSERT(a < numGenerators);
    POLY_ASSERT(b < numGenerators);

    vector<EdgeHash> meshEdges;
    for (set<unsigned>::const_iterator tetItr = tets.begin();
         tetItr != tets.end();
         ++tetItr) {
      i = *tetItr;
      POLY_ASSERT(i < out.numberoftetrahedra);
      POLY_ASSERT(tet2id.find(i) != tet2id.end());
      ii = tet2id[i];

      // Look for edges with adjacent tets.
      findOtherTetIndices(&out.tetrahedronlist[4*i], a, b, c, d);
      abc = hashFacet(a, b, c);
      abd = hashFacet(a, b, d);

      // Is abc a surface facet?
      if (facet2tets[abc].size() == 1) {
        POLY_ASSERT(facet2tets[abc][0] == i);
        POLY_ASSERT(facet2id.find(abc) != facet2id.end());
        jj = facet2id[abc];
        POLY_ASSERT(jj != ii);
        meshEdges.push_back(internal::hashEdge(ii, jj));
      } else {
        POLY_ASSERT((facet2tets[abc].size() == 2 and facet2tets[abc][0] == i) or facet2tets[abc][1] == i);
        k = (facet2tets[abc][0] == i ? facet2tets[abc][1] : facet2tets[abc][0]);
        jj = tet2id[k];
        if (jj != ii) meshEdges.push_back(internal::hashEdge(ii, jj));
      }

      // Is abd a surface facet?
      if (facet2tets[abd].size() == 1) {
        POLY_ASSERT(facet2tets[abd][0] == i);
        POLY_ASSERT(facet2id.find(abd) != facet2id.end());
        jj = facet2id[abd];
        POLY_ASSERT(jj != ii);
        meshEdges.push_back(internal::hashEdge(ii, jj));
      } else {
        POLY_ASSERT((facet2tets[abd].size() == 2) and (facet2tets[abd][0] == i or facet2tets[abd][1] == i));
        k = (facet2tets[abd][0] == i ? facet2tets[abd][1] : facet2tets[abd][0]);
        jj = tet2id[k];
        if (jj != ii) meshEdges.push_back(internal::hashEdge(ii, jj));
      }
    }

    // Arrange the edges in the correctly sorted and sign oriented order
    // to construct our face.
    sort(meshEdges.begin(), meshEdges.end());
    meshEdges.erase(unique(meshEdges.begin(), meshEdges.end()), meshEdges.end());
    if (meshEdges.size() > 1) {
      vector<int> edgeOrder;
      const bool infEdge = internal::computeSortedFaceEdges(meshEdges, edgeOrder);
      if (meshEdges.size() > 2) {

        // Add the edges and face to the quantized mesh.
        vector<int> face;
        for (vector<int>::const_iterator itr = edgeOrder.begin();
             itr != edgeOrder.end();
             ++itr) {
          const bool flip = (*itr < 0);
          k = (flip ? ~(*itr) : *itr);
          iedge = qmesh.addNewEdge(meshEdges[k]);
          face.push_back(flip ? ~iedge : iedge);
        }
        iface = qmesh.addNewFace(face);
        POLY_ASSERT(iface == qmesh.faces.size() - 1);

        // Add the face to its cells.
        e0 = qmesh.edgePosition(meshEdges[internal::positiveID(edgeOrder[0])]);
        e1 = qmesh.edgePosition(meshEdges[internal::positiveID(edgeOrder[1])]);
        e2 = qmesh.edgePosition(meshEdges[internal::positiveID(edgeOrder[2])]);
        vol = geometry::tetrahedralVolume6(&qmesh.generators[3*a], &e2.x, &e1.x, &e0.x);
        POLY_ASSERT(vol != 0.0);
        if (vol > 0.0) {
          qmesh.cells[a].push_back(iface);
          qmesh.cells[b].push_back(~iface);
        } else {
          qmesh.cells[a].push_back(~iface);
          qmesh.cells[b].push_back(iface);
        }

        // Did we create a new infEdge?  If so we know it was the second one in the ordered list.
        if (infEdge) {
          j = internal::positiveID(edgeOrder[1]);
          k = qmesh.edge2id[meshEdges[j]];
          qmesh.infEdges.push_back(k);
          cellInfEdges[a].push_back(meshEdges[j]);
          cellInfEdges[b].push_back(meshEdges[j]);
        }
      }
    }
  }

  // Build any infFaces we need.
  // For now we assume there is at most one infFace per cell, which is not true for all
  // degenerate cases!  Fix at some point...
  for (i = 0; i != numGenerators; ++i) {
    if (cellInfEdges[i].size() > 2) {
      vector<int> edgeOrder;
      internal::computeSortedFaceEdges(cellInfEdges[i], edgeOrder);

      // Check if we need to reverse the face node order.
      e0 = qmesh.edgePosition(cellInfEdges[i][internal::positiveID(edgeOrder[0])]);
      e1 = qmesh.edgePosition(cellInfEdges[i][internal::positiveID(edgeOrder[1])]);
      e2 = qmesh.edgePosition(cellInfEdges[i][internal::positiveID(edgeOrder[2])]);
      vol = geometry::tetrahedralVolume6(&qmesh.generators[3*a], &e2.x, &e1.x, &e0.x);
      POLY_ASSERT(vol != 0.0);
      if (vol < 0.0) {
        reverse(edgeOrder.begin(), edgeOrder.end());
        for (j = 0; j != edgeOrder.size(); ++j) edgeOrder[j] = ~edgeOrder[j];
      }
      iface = qmesh.faces.size();
      qmesh.faces.push_back(vector<int>());
      for (vector<int>::const_iterator itr = edgeOrder.begin();
           itr != edgeOrder.end();
           ++itr) {
        const bool flip = (*itr < 0);
        k = (flip ? ~(*itr) : *itr);
        iedge = qmesh.addNewEdge(cellInfEdges[i][k]);
        qmesh.faces[iface].push_back(flip ? ~iedge : iedge);
      }
      qmesh.cells[i].push_back(iface);
      qmesh.infFaces.push_back(iface);
    }
  }

  // Post-conditions.
  qmesh.assertValid();
}



//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
int64_t TetgenTessellator::coordMax = (1LL << 34);
double TetgenTessellator::mDegeneracy = 1.0/TetgenTessellator::coordMax;

}
