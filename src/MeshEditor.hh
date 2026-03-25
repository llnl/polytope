#ifndef POLYTOPE_MESHEDITOR_HH
#define POLYTOPE_MESHEDITOR_HH

#include <vector>
#include <iostream>
#include <algorithm>

#include "polytope.hh"
#include "Tessellation.hh"

namespace polytope {

// Typedef to define an edges in Polytope:
// Polytope doesn't store edge information separately. Define an EdgeType by the
// pair of (faceIndex, nodeNumber) where faceIndex is the positive mesh index
// for a face, and nodeNumber is the n-th node in the face array. The edge is
// then the segment connecting node n and n+1 (modulo numNodes) on that face.
// For 2D tessellation, nodeNumber is always 0, as faces and edges are equivalent
typedef std::pair<int, int> EdgeType;


// Comparison operator for Polytope edges
struct edgeCompare {
  bool operator()(const EdgeType e1, const EdgeType e2) const {
    return (e1.first <  e2.first                           ? true :
	    e1.first == e2.first and e1.second < e2.second ? true :
	    false);
  }
};


template<int Dimension, typename RealType>
class MeshEditor {

public:
  //-----------------------===== Public Interface =====-----------------------//

  // Constructors, destructors.
  MeshEditor(Tessellation<Dimension, RealType>& mesh);
  
  ~MeshEditor();

  // Deletes the mesh elements for the specified indices and recomputes the 
  // resulting mesh topology 
  void deleteCells(const std::vector<unsigned>& cellsToDelete);
  void deleteFaces(const std::vector<unsigned>& facesToDelete);
  void deleteNodes(const std::vector<unsigned>& nodesToDelete);

  // Clean small edges in the mesh based on some tolerance
  void cleanEdges(const RealType edgeTol);

  // Static variables
  static unsigned minEdgesPerFace;

private:
  //-----------------------===== Private Interface =====----------------------//

  // Compute face and node masks to clean small edges
  bool flagEdgesForCleaning(const RealType edgeTol,
                            std::vector<unsigned>& cellMap,
                            std::vector<unsigned>& faceMap,
                            std::vector<unsigned>& nodeMap);
  
  // Does the work of recomputing mesh based on cell/face/node masks
  void cleanMesh(std::vector<unsigned>& cellMap,
                 std::vector<unsigned>& faceMap,
                 std::vector<unsigned>& nodeMap);

  // First computes the new->old element maps then passes to routine above
  void cleanMesh();
  
  // Used by the delete[Element] routines to build up the necessary element
  // masks used by cleanMesh()
  void computeMasks();

  // Computes the other two element masks, given one
  void computeFaceAndNodeMasks();
  void computeCellAndNodeMasks();
  void computeCellAndFaceMasks();


  //-----------------------===== Class Members =====--------------------------//

  // Hold a reference to the mesh
  Tessellation<Dimension, RealType>& mMesh;

  // The element masks
  std::vector<unsigned> mCellMask, mFaceMask, mNodeMask;

  // Disallowed
  MeshEditor();
  MeshEditor(const MeshEditor&);
};

}// end polytope namespace

#endif
