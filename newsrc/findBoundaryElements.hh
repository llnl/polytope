//------------------------------------------------------------------------------
// Identify the boundary nodes and faces in the given tessellation.
// This method uses a topological method -- if you walk all the faces of the
// cells in a tessellation, you'll hit the boundary faces once and the interior
// faces twice.
//------------------------------------------------------------------------------
#ifndef __polytope_findBoundaryElements__
#define __polytope_findBoundaryElements__

#include <vector>
#include <algorithm>

#include "Tessellation.hh"
#include "polytope_internal.hh"


namespace polytope {
template<int Dimension, typename RealType>
void
findBoundaryElements(const Tessellation<Dimension, RealType>& mesh,
                     std::vector<unsigned>& boundaryFaces,
                     std::vector<unsigned>& boundaryNodes) {
  boundaryFaces.clear();
  boundaryNodes.clear();
  for (unsigned iface = 0; iface < mesh.faces.size(); ++iface) {
    POLY_ASSERT2(mesh.faceCells[iface].size() == 1 or
                 mesh.faceCells[iface].size() == 2,
                 mesh.faceCells[iface].size());
    if (mesh.faceCells[iface].size() == 1) {
      boundaryFaces.push_back(iface);
      std::copy(mesh.faces[iface].begin(), mesh.faces[iface].end(), std::back_inserter(boundaryNodes));
    }
  }
  std::sort(boundaryNodes.begin(), boundaryNodes.end());
  boundaryNodes.erase(unique(boundaryNodes.begin(), boundaryNodes.end()), boundaryNodes.end());

  POLY_ASSERT(boundaryNodes.size() > 0);
  // Note: In 2D, boundary nodes ~ O(sqrt(N)). In 3D, boundary nodes ~ O(N^(2/3)).
  // For small meshes, thin domains, or surface tessellations, boundary nodes can be
  // a large fraction of total nodes. We just verify there are some boundary nodes.
  POLY_ASSERT(boundaryNodes.size() <= mesh.nodes.size()/Dimension);
}

}

#else

namespace polytope {
template<int Dimension, typename RealType>
void
findBoundaryElements(const Tessellation<Dimension, RealType>& mesh,
                     std::vector<int>& boundaryFaces,
                     std::vector<int>& boundaryNodes);
}

#endif

