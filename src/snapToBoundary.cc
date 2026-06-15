//------------------------------------------------------------------------------
// Take a tessellation and a PLC.  Make sure all points within the given
// tolerance are on the boundary, and that all PLC boundary points are accounted
// for.
//------------------------------------------------------------------------------
#include <set>

#include "snapToBoundary.hh"
#include "Point.hh"
#include "KeyTraits.hh"
#include "nearestPoint.hh"
#include "polytope_tessellator_utilities.hh"

using namespace std;

namespace polytope {

//------------------------------------------------------------------------------
// 2D
//------------------------------------------------------------------------------
template<typename RealType>
void snapToBoundary(Tessellation<2, RealType>& mesh,
                    const vector<RealType>& points,
                    const PLC<2>& geometry,
                    const RealType degeneracy) {

  typedef Point2<RealType> RealPoint;

  // Find the bounding box, and our tolerance for points being on the boundary.
  RealPoint xmin, xmax;
  RealType length;
  geometry::computeBoundingBox<2, RealType>(&points[0], points.size(), true, &xmin.x, &xmax.x);
  length = std::max(xmax.x - xmin.x, xmax.y - xmin.y);
  const RealType tol = 4.0*degeneracy*length;

  // Copy the boundary nodes to a std::set.
  std::set<unsigned> boundaryNodes(mesh.boundaryNodes.begin(), mesh.boundaryNodes.end());

  // Snap a subset of the boundary nodes to the PLC boundary locations
  set<unsigned> plcNodes;
  set<unsigned> indices;
  for (unsigned ii = 0; ii < points.size()/2; ++ii) indices.insert(ii);
  for (set<unsigned>::iterator iitr = indices.begin();  
       iitr != indices.end();) {
    unsigned ipoint = *iitr;
    unsigned inode;
    const RealPoint rp = RealPoint(points[2*ipoint], points[2*ipoint+1]);
    RealType dist = std::numeric_limits<RealType>::max();
    set<unsigned>::iterator nodeItr = boundaryNodes.begin();
    while (nodeItr != boundaryNodes.end() and dist > tol) {
      inode = *nodeItr;
      dist = geometry::distance<2, RealType>(&mesh.nodes[2*(*nodeItr)], &rp.x);
      ++nodeItr;
    }
    if (dist < tol) {
      iitr++;
      plcNodes.insert(inode);
      mesh.nodes[2*inode  ] = rp.x;
      mesh.nodes[2*inode+1] = rp.y;
      boundaryNodes.erase(inode);
      indices.erase(ipoint);
    } else {
      ++iitr;
    }
  }
 
  // Check that all of the PLC points are accounted for
  POLY_BEGIN_CONTRACT_SCOPE;
  {
    if (not indices.empty()) {
      for (set<unsigned>::iterator itr = indices.begin();
           itr != indices.end();
           ++itr) {
        RealType dist = std::numeric_limits<RealType>::max();
        unsigned inode;
        for (unsigned ii = 0; ii < mesh.nodes.size()/2; ++ii) {
          RealType dnew = geometry::distance<2,RealType>(&points[2*(*itr)], &mesh.nodes[2*ii]);
          if (dnew < dist) {
            inode = ii;
            dist = dnew;
          }
        }
        cerr << "PLC point " << *itr << endl
             << "  location = (" << points[2*(*itr)  ] << "," 
             << points[2*(*itr)+1] << ")" << endl
             << "  nearNode = " << inode << endl
             << "  node loc = (" << mesh.nodes[2*inode] << "," 
             << mesh.nodes[2*inode+1] << ")" << endl
             << "  distance = " << dist << "      tol = " << tol << endl
             << "  boundaryNode? " 
             << (boundaryNodes.find(inode) != boundaryNodes.end() ? "YES" : "NO") 
             << endl;
      }
    }
    POLY_ASSERT(indices.empty());
    POLY_ASSERT2(plcNodes.size() == points.size()/2,
                 plcNodes.size() << " != " << points.size()/2);
  }
  POLY_END_CONTRACT_SCOPE;

  // Snap all boundary nodes to the nearest floating point along boundary
  for (set<unsigned>::iterator nodeItr = boundaryNodes.begin();
       nodeItr != boundaryNodes.end();
       ++nodeItr) {
    RealPoint result;
    RealType dist = nearestPoint(&mesh.nodes[2*(*nodeItr)],
                                 points.size()/2,
                                 &points[0],
                                 geometry,
                                 &result.x);
    // POLY_ASSERT2(dist < tol,
    //              "Possible internal boundary node "
    //              << (*nodeItr) << " at ("
    //              << mesh.nodes[2*(*nodeItr)  ] << ","
    //              << mesh.nodes[2*(*nodeItr)+1] << ")"
    //              << " wants to move distance " << dist << " to "
    //              << "(" << result[0] << "," << result[1]
    //              << ") over tolerance " << tol);
    mesh.nodes[2*(*nodeItr)  ] = result.x;
    mesh.nodes[2*(*nodeItr)+1] = result.y;
  }
}

//------------------------------------------------------------------------------
// 3D
//------------------------------------------------------------------------------
template<typename RealType>
void snapToBoundary(Tessellation<3, RealType>& mesh,
                    const vector<RealType>& points,
                    const PLC<3>& geometry,
                    const RealType degeneracy) {
  (void) degeneracy;
  (void) geometry;
  (void) points;
  (void) mesh;
  POLY_ASSERT2(false, "Implement me!");
}

//------------------------------------------------------------------------------
// Instantiation.
//------------------------------------------------------------------------------
template void snapToBoundary<double>(Tessellation<2, double>& mesh,
                                     const vector<double>& points,
                                     const PLC<2>& geometry,
                                     const double degeneracy);
template void snapToBoundary<double>(Tessellation<3, double>& mesh,
                                     const vector<double>& points,
                                     const PLC<3>& geometry,
                                     const double degeneracy);

}
