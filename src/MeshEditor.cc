//------------------------------------------------------------------------
// MeshEditor
//------------------------------------------------------------------------
#include <limits>
#include <numeric>

#include "polytope.hh"
#include "MeshEditor.hh"
#include "polytope_geometric_utilities.hh"

namespace polytope {

using namespace std;

//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
MeshEditor<Dimension, RealType>::
MeshEditor(Tessellation<Dimension, RealType>& mesh):
  mMesh(mesh){
  minEdgesPerFace = (Dimension == 2) ? 1 : 3;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
MeshEditor<Dimension, RealType>::
~MeshEditor(){
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
MeshEditor<Dimension, RealType>::
deleteCells(const std::vector<unsigned>& cellsToDelete) {
  // Pre-conditions
  POLY_ASSERT2(!cellsToDelete.empty(), "No cells specified by deletion");
  const unsigned ncells = mMesh.cells.size();
  mCellMask = std::vector<unsigned>(ncells, 1);
  for (std::vector<unsigned>::const_iterator itr = cellsToDelete.begin();
       itr != cellsToDelete.end(); ++itr) {
    POLY_ASSERT(*itr < ncells);
    mCellMask[*itr] = 0;
  }
  computeMasks();
  cleanMesh();
  // Post-conditions
  POLY_ASSERT(mCellMask.empty() and mFaceMask.empty() and mNodeMask.empty());
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
MeshEditor<Dimension, RealType>::
deleteFaces(const std::vector<unsigned>& facesToDelete) {
  POLY_ASSERT2(false, "This routine is currently not implemented");  
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
MeshEditor<Dimension, RealType>::
deleteNodes(const std::vector<unsigned>& nodesToDelete) {
  POLY_ASSERT2(false, "This routine is currently not implemented");
}
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
MeshEditor<Dimension, RealType>::
cleanEdges(const RealType edgeTol) {
  // Pre-conditions
  POLY_ASSERT2(edgeTol > 0.0, "Specify a positive (non-zero) edge tolerance!");

  // flag edges smaller than the given tolerance in terms of face and node masks
  std::vector<unsigned> cellMap, faceMap, nodeMap;
  bool edgesClean = flagEdgesForCleaning(edgeTol, cellMap, faceMap, nodeMap);
  if (!edgesClean) cleanMesh(cellMap, faceMap, nodeMap);

  // Post-conditions
  POLY_ASSERT(mCellMask.empty() and mFaceMask.empty() and mNodeMask.empty());
}
//------------------------------------------------------------------------------



//-----------------------===== Private Methods =====------------------------//



//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
bool
MeshEditor<Dimension, RealType>::
flagEdgesForCleaning(const RealType edgeTol,
                     vector<unsigned>& cellMap,
                     vector<unsigned>& faceMap,
                     vector<unsigned>& nodeMap) {

  // Pre-conditions.
  POLY_ASSERT2(edgeTol > 0.0, "Specify a positive (non-zero) edge tolerance!");

  bool edgesClean = true;
  const unsigned ncells0 = mMesh.cells.size();
  const unsigned nfaces0 = mMesh.faces.size();
  const unsigned nnodes0 = mMesh.nodes.size()/Dimension;

  // Map edges to a mesh ID, compute maximum edge lengths per cell, and map cells
  // to edges
  // //EdgeMap edgeToID;
  map<EdgeType, int> edgeToIndex;
  map<int, EdgeType> indexToEdge;
  map<unsigned, vector<unsigned> > cellToEdges, borderNodesToBorderFaces;
  set<unsigned> borderNodes;
  vector<RealType> edgeLength;
  RealType length;
  unsigned inode0, inode1;
  // //int edgeCount = 0;
  for (unsigned icell = 0; icell != ncells0; ++icell) {
    vector<unsigned> cellEdgeIDs;
    for (vector<int>::const_iterator itr = mMesh.cells[icell].begin();
	 itr != mMesh.cells[icell].end(); ++itr) {
      const unsigned iface = (*itr < 0) ? ~(*itr) : *itr;
      POLY_ASSERT(iface < mMesh.faces.size());
      const unsigned nfaceNodes = mMesh.faces[iface].size();
      POLY_ASSERT(nfaceNodes >= 2);
      const unsigned nfaceCells = mMesh.faceCells[iface].size();
      POLY_ASSERT(nfaceCells == 1 or nfaceCells == 2);
      if (nfaceCells == 1) {
        for (std::vector<unsigned>::iterator nodeItr = mMesh.faces[iface].begin();
             nodeItr != mMesh.faces[iface].end(); ++nodeItr) {
          borderNodes.insert(*nodeItr);
          borderNodesToBorderFaces[*nodeItr].push_back(iface);
        }
      }
      const unsigned maxNodeIndex = (nfaceNodes == 2) ? 1 : nfaceNodes;
      for (unsigned inode = 0; inode != maxNodeIndex; ++inode) {
        inode0 = mMesh.faces[iface][inode];
        inode1 = mMesh.faces[iface][(inode+1) % nfaceNodes];
        POLY_ASSERT(inode0 != inode1);
        POLY_ASSERT(inode0 < mMesh.nodes.size()/Dimension and
                    inode1 < mMesh.nodes.size()/Dimension);
        const EdgeType edge = internal::hashEdge(inode0, inode1);
        const unsigned nedges = edgeToIndex.size();
        // const int iedge = internal::addKeyToMap<EdgeType,internal::pairCompare<unsigned, unsigned> >(edge, edgeToIndex);
        const int iedge = internal::addKeyToMap(edge, edgeToIndex);
        if (iedge == nedges) {
           indexToEdge[iedge] = edge;
           length = geometry::distance<Dimension, RealType>
             (&mMesh.nodes[Dimension*inode0], &mMesh.nodes[Dimension*inode1]);
           edgeLength.push_back(length);
        }
        cellEdgeIDs.push_back(iedge);

        // // typename EdgeMap::left_const_iterator lItr = edgeToID.left.find(edge);
	// // if (lItr == edgeToID.left.end()) {
	// //   edgeToID.insert(typename EdgeMap::value_type(edge, edgeCount));
	// //   length = geometry::distance<Dimension, RealType>
	// //     (&mMesh.nodes[Dimension*inode0], &mMesh.nodes[Dimension*inode1]);
	// //   edgeLength.push_back(length);
	// //   cellEdgeIDs.push_back(edgeCount);
	// //   ++edgeCount;
	// // } else {
	// //   cellEdgeIDs.push_back(lItr->second);
	// // }
      }
    }
    cellToEdges[icell] = cellEdgeIDs;
  }
  POLY_ASSERT(edgeToIndex.size() == edgeLength.size());
  POLY_ASSERT(cellToEdges.size() == ncells0);
  POLY_ASSERT(edgeToIndex.size() == indexToEdge.size());
  const unsigned nedges0 = edgeToIndex.size();

  // Create a simple lookup for the border nodes. Also create a lookup
  // for corner nodes: nodes corresponding to 
  set<unsigned> cornerNodes;
  // vector<unsigned> isBorderNode(nnodes0, 0), isCornerNode(nnodes0, 0);
  for (set<unsigned>::const_iterator nodeItr = borderNodes.begin();
       nodeItr != borderNodes.end(); 
       ++nodeItr) {
    POLY_ASSERT(*nodeItr < nnodes0);
    // isBorderNode[*nodeItr] = 1;
    vector<unsigned> borderFaces = borderNodesToBorderFaces[*nodeItr];
    POLY_ASSERT(borderFaces.size() == 2);
    vector<unsigned> otherNodes;
    for (vector<unsigned>::const_iterator faceItr = borderFaces.begin();
         faceItr != borderFaces.end();
         ++faceItr) {
      const unsigned iface = *faceItr;
      POLY_ASSERT(iface < mMesh.faces.size());
      POLY_ASSERT(mMesh.faceCells[iface].size() == 1);
      POLY_ASSERT(mMesh.faces[iface].size() == 2);
      const unsigned otherNode = (mMesh.faces[iface][0] == *nodeItr) ?
         mMesh.faces[iface][1] : mMesh.faces[iface][0];
      otherNodes.push_back(otherNode);
    }
    POLY_ASSERT(otherNodes.size() == 2);

    const bool collinear =
      geometry::collinear<Dimension, RealType>(&mMesh.nodes[2*(*nodeItr)],
                                               &mMesh.nodes[2*otherNodes[0]],
                                               &mMesh.nodes[2*otherNodes[1]],
                                               1.0e-8);
    if (not collinear) cornerNodes.insert(*nodeItr);// isCornerNode[*nodeItr] = 1;
  }
  
  // Compute the maximum edge length for the cells around an edge
  vector<RealType> maxCellEdgeLength(nedges0, 0.0);
  for (unsigned icell = 0; icell != ncells0; ++icell) {
    RealType maxCellEdge = 0.0;
    vector<unsigned> edgeIDs = cellToEdges[icell];
    for (vector<unsigned>::const_iterator itr = edgeIDs.begin();
	 itr != edgeIDs.end();
	 ++itr) maxCellEdge = max(maxCellEdge, edgeLength[*itr]);
    for (vector<unsigned>::const_iterator itr = edgeIDs.begin();
	 itr != edgeIDs.end();
	 ++itr) maxCellEdgeLength[*itr] = max(maxCellEdgeLength[*itr], maxCellEdge);
  }
  
  // Flag the edges we want to remove.
  vector<unsigned> edgeMask(nedges0, 1);
  vector<unsigned> nodeCollapse(nnodes0);
  EdgeType edge;
  mNodeMask = vector<unsigned>(nnodes0, 1);
  for (unsigned i = 0; i != nnodes0; ++i) nodeCollapse[i] = i;
  for (unsigned iedge = 0; iedge != nedges0; ++iedge) {
    POLY_ASSERT(indexToEdge.find(iedge) != indexToEdge.end());
    edge = indexToEdge[iedge];

    // // typename EdgeMap::right_const_iterator itr = edgeToID.right.find(iedge);
    // // POLY_ASSERT(itr != edgeToID.right.end());
    // // edge = itr->second;

    inode0 = edge.first;
    inode1 = edge.second;
    const bool cornerEdge = (cornerNodes.find(inode0) != cornerNodes.end() or
                             cornerNodes.find(inode1) != cornerNodes.end());
                            //(isCornerNode[inode0] == 0 or isCornerNode[inode1] == 0));
    const bool cleanTest = (edgeLength[iedge] < edgeTol*maxCellEdgeLength[iedge] and
                            mNodeMask[inode0] == 1 and 
                            mNodeMask[inode1] == 1 and
                            not cornerEdge);
    if (cleanTest) {
      edgesClean = false;
      edgeMask[iedge] = 0;
      const bool node1KeepTest = //isCornerNode[inode1] == 1 or 
         (cornerNodes.find(inode1) != cornerNodes.end()) or
         (borderNodes.find(inode1) != borderNodes.end() and
          borderNodes.find(inode0) == borderNodes.end());
         //(isBorderNode[inode1] == 1 and isBorderNode[inode0] == 0);
      unsigned keepNode   = node1KeepTest ? inode1 : inode0;
      unsigned deleteNode = node1KeepTest ? inode0 : inode1;
      mNodeMask[keepNode  ] = 2;
      mNodeMask[deleteNode] = 0;
      nodeCollapse[deleteNode] = keepNode;
    }
  }
  replace_if(mNodeMask.begin(), mNodeMask.end(), bind2nd(equal_to<unsigned>(), 2), 1);

  
#ifdef POLYTOPE_ENABLE_MPI
  // Parallel configuration.
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // TODO: Put in an all-reduced switch based on the edgesClean flag. After all, if
  //       all processors have clean mesh, why bother communicating?

  // Construct a list of shared edges, sorted the same on each proc. Construct and 
  // communicate the masks and make the local edge mask consistent across domains.
  // Of course, none of this is necessary if we're running in serial, so...
  if (numProcs > 1) {

    // // Create reverse loopup from edge to collection of neighboring cells
    // map<EdgeType, vector<unsigned> > edgeToCells;
    // for (map<unsigned, vector<unsigned> >::const_iterator cellItr = cellToEdges.begin();
    //      cellItr != cellToEdges.end(); ++cellItr) {
    //   vector<unsigned> 
    // }

    // Construct shared edge list
    vector<vector<EdgeType> > sharedEdges   (mMesh.neighborDomains.size());
    vector<vector<unsigned> > sharedEdgeMask(mMesh.neighborDomains.size());
    vector<vector<unsigned> > sharedEdgeReverse(nedges0);
    vector<vector<unsigned> > sharedNodeReverse(nnodes0);
    vector<map<unsigned,unsigned> > domainNodeMaps(mMesh.neighborDomains.size());
    unsigned ineighbor = 0;
    for (vector<unsigned>::const_iterator otherItr = mMesh.neighborDomains.begin();
         otherItr != mMesh.neighborDomains.end(); ++otherItr, ++ineighbor) {
      
      // Create reverse lookup from nodeID to position in sharedNodes list
      map<unsigned, unsigned> nodeIndexToSharePosition;
      for (vector<unsigned>::const_iterator nodeItr = mMesh.sharedNodes[ineighbor].begin();
           nodeItr != mMesh.sharedNodes[ineighbor].end(); ++nodeItr) {
        POLY_ASSERT(nodeIndexToSharePosition.find(*nodeItr) == nodeIndexToSharePosition.end());
        unsigned result = nodeIndexToSharePosition.size();
        nodeIndexToSharePosition[*nodeItr] = result;
        sharedNodeReverse[*nodeItr].push_back(ineighbor);
      }
      POLY_ASSERT(nodeIndexToSharePosition.size() == mMesh.sharedNodes[ineighbor].size());
      domainNodeMaps[ineighbor] = nodeIndexToSharePosition;

      // Consider all edges on this domain. The shared ones consist of two shared nodes.
      for (typename map<EdgeType, int>::const_iterator edgeItr = edgeToIndex.begin();
           edgeItr != edgeToIndex.end();
           ++edgeItr) {

      // // for (typename EdgeMap::left_const_iterator edgeItr = edgeToID.left.begin();
      // //      edgeItr != edgeToID.left.end(); ++edgeItr) {

        const EdgeType edge  = edgeItr->first;
        const unsigned iedge = edgeItr->second;
        inode0 = edge.first;
        inode1 = edge.second;
        POLY_ASSERT(inode0 < inode1);
        POLY_ASSERT(inode0 < mMesh.nodes.size()/Dimension and
                    inode1 < mMesh.nodes.size()/Dimension);
        map<unsigned,unsigned>::const_iterator it0 = nodeIndexToSharePosition.find(inode0);
        map<unsigned,unsigned>::const_iterator it1 = nodeIndexToSharePosition.find(inode1);
        if (it0 != nodeIndexToSharePosition.end() and 
            it1 != nodeIndexToSharePosition.end()) {
          const EdgeType sharedEdge = internal::hashEdge(it0->second, it1->second);
          sharedEdges[ineighbor].push_back(sharedEdge);
          sharedEdgeReverse[iedge].push_back(ineighbor);
        }
      }

      // Sort the sharedEdges list
      sort(sharedEdges[ineighbor].begin(), sharedEdges[ineighbor].end(),
           internal::pairCompare<unsigned, unsigned>);        

      // // Blago!
      // cerr << rank << " <--> " << mMesh.neighborDomains[ineighbor] << " :" << endl;
      // cerr << "   Nodes: ";
      // copy(mMesh.sharedNodes[ineighbor].begin(), mMesh.sharedNodes[ineighbor].end(), 
      //      ostream_iterator<unsigned>(cerr, " "));
      // cerr << endl << "   Edges:";
      // for (unsigned ii = 0; ii != sharedEdges[ineighbor].size(); ++ii) {
      //    cerr << " (" << sharedEdges[ineighbor][ii].first 
      //         << ","  << sharedEdges[ineighbor][ii].second << ")";
      // }
      // cerr << endl;
      // // Blago!


      // Now make the corresponding mask for each sharedEdge list
      sharedEdgeMask[ineighbor].resize(sharedEdges[ineighbor].size());
      unsigned edgeCount = 0;
      for (vector<EdgeType>::const_iterator sEdgeItr = sharedEdges[ineighbor].begin();
           sEdgeItr != sharedEdges[ineighbor].end(); ++sEdgeItr, ++edgeCount) {
        const unsigned ind0 = sEdgeItr->first;
        const unsigned ind1 = sEdgeItr->second;
        POLY_ASSERT(ind0 < mMesh.sharedNodes[ineighbor].size() and
                    ind1 < mMesh.sharedNodes[ineighbor].size());
        inode0 = mMesh.sharedNodes[ineighbor][ind0];
        inode1 = mMesh.sharedNodes[ineighbor][ind1];
        POLY_ASSERT(inode0 < mMesh.nodes.size()/Dimension and
                    inode1 < mMesh.nodes.size()/Dimension);
        const EdgeType edge = internal::hashEdge(inode0, inode1);

        typename map<EdgeType, int>::const_iterator itr = edgeToIndex.find(edge);
        POLY_ASSERT(itr != edgeToIndex.end());
        
        // // typename EdgeMap::left_const_iterator itr = edgeToID.left.find(edge);
        // // POLY_ASSERT(itr != edgeToID.left.end());

        const unsigned iedge = itr->second;
        POLY_ASSERT(iedge < nedges0);
        sharedEdgeMask[ineighbor][edgeCount] = edgeMask[iedge];
      }
      POLY_ASSERT(edgeCount == sharedEdges[ineighbor].size());
    }
    
    // Send the sharedEdgeMask to each neighbor asynchronously
    vector<char> localBuffer;
    vector<MPI_Request> sendRequests;
    sendRequests.reserve(2*mMesh.neighborDomains.size());
    ineighbor = 0;
    for (vector<unsigned>::const_iterator otherItr = mMesh.neighborDomains.begin();
         otherItr != mMesh.neighborDomains.end(); ++otherItr, ++ineighbor) {
      serialize(sharedEdgeMask[ineighbor], localBuffer);
      unsigned localBufferSize = localBuffer.size();
      const unsigned otherProc = *otherItr;
      sendRequests.push_back(MPI_Request());
      MPI_Isend(&localBufferSize, 1, MPI_UNSIGNED, otherProc, 1, MPI_COMM_WORLD, &sendRequests.back());
      if (localBufferSize > 0) {
        sendRequests.push_back(MPI_Request());
        MPI_Isend(&localBuffer.front(), localBufferSize, MPI_CHAR, otherProc, 2, MPI_COMM_WORLD, &sendRequests.back());
      }
    }
    POLY_ASSERT(sendRequests.size() <= 2*mMesh.neighborDomains.size());

    // Get the masks from the neighbors
    vector<vector<unsigned> > neighborEdgeMasks(mMesh.neighborDomains.size());
    ineighbor = 0;
    for (vector<unsigned>::const_iterator otherItr = mMesh.neighborDomains.begin();
         otherItr != mMesh.neighborDomains.end(); ++otherItr, ++ineighbor) {
      const unsigned otherProc = *otherItr;
      unsigned bufSize;
      MPI_Status recvStatus1, recvStatus2;
      MPI_Recv(&bufSize, 1, MPI_UNSIGNED, otherProc, 1, MPI_COMM_WORLD, &recvStatus1);
      if (bufSize > 0) {
        vector<char> buffer(bufSize, '\0');
        MPI_Recv(&buffer.front(), bufSize, MPI_CHAR, otherProc, 2, MPI_COMM_WORLD, &recvStatus2);
        vector<char>::const_iterator itr = buffer.begin();
        vector<unsigned> otherEdgeMask;
        deserialize(otherEdgeMask, itr, buffer.end());
        POLY_ASSERT(itr == buffer.end());
        neighborEdgeMasks[ineighbor] = otherEdgeMask;
      }
    }
    
    // Make sure all the sends are complete
    vector<MPI_Status> sendStatus(sendRequests.size());
    MPI_Waitall(sendRequests.size(), &sendRequests.front(), &sendStatus.front());
    
    // Make the edge masks consistent across domain boundaries. If one domain 
    // says delete (mask=0) and one says keep (mask=1), we always keep.
    for (unsigned i = 0; i != mMesh.neighborDomains.size(); ++i) {
      unsigned edgeCount = 0;
      for (vector<EdgeType>::const_iterator sEdgeItr = sharedEdges[i].begin();
           sEdgeItr != sharedEdges[i].end(); ++sEdgeItr, ++edgeCount) {
        const unsigned ind0 = sEdgeItr->first;
        const unsigned ind1 = sEdgeItr->second;
        POLY_ASSERT(ind0 < mMesh.sharedNodes[i].size() and
                    ind1 < mMesh.sharedNodes[i].size());
        inode0 = mMesh.sharedNodes[i][ind0];
        inode1 = mMesh.sharedNodes[i][ind1];
        POLY_ASSERT(inode0 < mMesh.nodes.size()/Dimension and
                    inode1 < mMesh.nodes.size()/Dimension);
        const EdgeType edge = internal::hashEdge(inode0, inode1);
        typename map<EdgeType, int>::const_iterator itr = edgeToIndex.find(edge);
        POLY_ASSERT(itr != edgeToIndex.end());

        // // typename EdgeMap::left_const_iterator itr = edgeToID.left.find(edge);
        // // POLY_ASSERT(itr != edgeToID.left.end());

        const unsigned iedge = itr->second;
        POLY_ASSERT(iedge < nedges0);
        //POLY_ASSERT(edgeCount < neighborEdgeMasks[i].size());
        edgeMask[iedge] = max(edgeMask[iedge], neighborEdgeMasks[i][edgeCount]);
      }
      POLY_ASSERT(edgeCount == sharedEdges[i].size());
    }

    // Ammend the edgesClean flag if we're not cleaning after all.
    if (std::accumulate(edgeMask.begin(), edgeMask.end(), 0) == nedges0)
       edgesClean = true;

    // We've invalidated the nodeCollapse and nodeMask lists. Make them consistent.
    mNodeMask = vector<unsigned>(nnodes0, 1);
    if (!edgesClean) {
      for (unsigned i = 0; i != nnodes0; ++i) nodeCollapse[i] = i;
      for (unsigned iedge = 0; iedge != nedges0; ++iedge) {
        typename map<int, EdgeType>::const_iterator itr = indexToEdge.find(iedge);
        POLY_ASSERT(itr != indexToEdge.end());

        // // typename EdgeMap::right_const_iterator itr = edgeToID.right.find(iedge);
        // // POLY_ASSERT(itr != edgeToID.right.end());

        if (edgeMask[iedge] == 0) {
          const EdgeType edge = itr->second;
          const bool isShared = (sharedEdgeReverse[iedge].size() == 0) ? false : true;
          inode0 = edge.first;
          inode1 = edge.second;
          const bool node0Shared = (sharedNodeReverse[inode0].size() == 0) ? false : true;
          const bool node1Shared = (sharedNodeReverse[inode1].size() == 0) ? false : true;

          if (isShared) {
            POLY_ASSERT2(sharedEdgeReverse[iedge].size() == 1, "I think this logic only"
                         << " holds for 2D meshes. An edge is a face and can be shared by"
                         << " at most two domains. In 3D, 3 or more domains can share an"
                         << " edge. It's not clear if the orderings of the shared node"
                         << " lists are consistent among 3 or more domains");
            const unsigned ineighbor = sharedEdgeReverse[iedge][0];
            const unsigned ind0 = domainNodeMaps[ineighbor][inode0];
            const unsigned ind1 = domainNodeMaps[ineighbor][inode1];
            if (ind1 < ind0) {inode0 = edge.second; inode1 = edge.first;}
          }

          else if (node1Shared and !node0Shared) {
             inode0 = edge.second; inode1 = edge.first;
          }

          mNodeMask[inode1] = 0;
          nodeCollapse[inode1] = inode0;

          // // Blago!
          // cerr << "Delete edge " << iedge << " having nodes "
          //      << "(" << edge.first << "," << edge.second << ")\n"
          //      << "    Is it shared? " << (isShared ? "Yes" : "NO") << "\n"
          //      << "    Delete " << inode1 << " and map its index to " << inode0 << endl;
          // // Blago!
        }
      }
    }

  } // Done with the multi-processor stuff
#endif
  
  // Check if there are any faces that need to be removed.
  // (Really only for 3D tessellations since edges=faces in 2D.)
  mFaceMask = std::vector<unsigned>(nfaces0, 1);
  if (!edgesClean) {
    unsigned numActiveEdges;
    for (unsigned iface = 0; iface != nfaces0; ++iface) {
      numActiveEdges = 0;
      const unsigned maxNodeIndex = (mMesh.faces[iface].size() == 2) ? 1 : mMesh.faces[iface].size();
      for (unsigned inode = 0; inode != maxNodeIndex; ++inode) {
        inode0 = mMesh.faces[iface][inode];
        inode1 = mMesh.faces[iface][(inode+1) % mMesh.faces[iface].size()];
        const EdgeType edge = internal::hashEdge(inode0, inode1);
        typename map<EdgeType, int>::const_iterator itr = edgeToIndex.find(edge);
        POLY_ASSERT(itr != edgeToIndex.end());
	numActiveEdges += edgeMask[itr->second];
      }

      // // Blago!
      // cerr << "Face " << iface << ": " 
      //      << "Num active edges = " << numActiveEdges << ": "
      //      << (numActiveEdges >= minEdgesPerFace ? "KEEP!" : "EXTERMINATE!") << endl;
      // // Blago!
     
      if( numActiveEdges < minEdgesPerFace) mFaceMask[iface] = 0;
    }
  }
  
  // Build up the node map from original indices to new indices
  unsigned nodeCount = 0;
  nodeMap.resize(nnodes0);
  for (unsigned i = 0; i != nnodes0; ++i) {
    nodeMap[i] = nodeCount;
    nodeCount += mNodeMask[i];
  }
  for (unsigned i = 0; i != nnodes0; ++i) {
    nodeMap[i] = nodeMap[nodeCollapse[i]];
  }
  POLY_ASSERT(nodeCount == std::accumulate(mNodeMask.begin(), mNodeMask.end(), 0));

  // Build up the face map from original indices to new indices
  faceMap.resize(nfaces0);
  unsigned faceCount = 0;
  for (unsigned iface = 0; iface != nfaces0; ++iface) {
    faceMap[iface] = faceCount;
    faceCount += mFaceMask[iface];
  }
  POLY_ASSERT(faceCount <= nfaces0);

  // // Blago!
  // cerr << "Delete faces: ";
  // for (unsigned i = 0; i != nfaces0; ++i) if(mFaceMask[i]==0) cerr << i << " ";
  // cerr << endl << "Delete nodes: ";
  // for (unsigned i = 0; i != nnodes0; ++i) if(mNodeMask[i]==0) cerr << i << " ";
  // cerr << endl;
  // // Blago!

  // Deleting edges (hopefully) does not result in a deleted cell
  mCellMask = std::vector<unsigned>(ncells0, 1);
  cellMap.resize(ncells0);
  for (unsigned i = 0; i != ncells0; ++i) cellMap[i] = i;

  //Post-conditions
  if (edgesClean) {
     POLY_ASSERT(std::accumulate(mNodeMask.begin(), mNodeMask.end(), 0) == nnodes0);
     POLY_ASSERT(std::accumulate(mFaceMask.begin(), mFaceMask.end(), 0) == nfaces0);
     POLY_ASSERT(std::accumulate(mCellMask.begin(), mCellMask.end(), 0) == ncells0);
     mNodeMask.clear();  mFaceMask.clear();  mCellMask.clear();
  } else {
    POLY_ASSERT(std::accumulate(mNodeMask.begin(), mNodeMask.end(), 0) < nnodes0);
  }
  
  return edgesClean;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
MeshEditor<Dimension, RealType>::
cleanMesh() {

  // Pre-conditions
  POLY_ASSERT(!mCellMask.empty() and !mFaceMask.empty() and !mNodeMask.empty());

  // Original mesh sizes
  const unsigned ncells0 = mMesh.cells.size();
  const unsigned nfaces0 = mMesh.faces.size();
  const unsigned nnodes0 = mMesh.nodes.size()/Dimension;

  // Determine the new cell, face, and node numberings.
  unsigned nnodes1 = 0, nfaces1 = 0, ncells1 = 0;
  std::vector<unsigned> nodeMap(nnodes0), faceMap(nfaces0), cellMap(ncells0);
  for (unsigned i = 0; i != nnodes0; ++i) {
    nodeMap[i] = nnodes1;
    if (mNodeMask[i] == 1) nnodes1++;
  }
  for (unsigned i = 0; i != nfaces0; ++i) {
    faceMap[i] = nfaces1;
    if (mFaceMask[i] == 1) nfaces1++;
  }
  for (unsigned i = 0; i != ncells0; ++i) {
    if (mCellMask[i] == 1) cellMap[i] = ncells1++;
  }

  // Post-conditions
  POLY_ASSERT(std::accumulate(mNodeMask.begin(), mNodeMask.end(), 0) == nnodes1);
  POLY_ASSERT(std::accumulate(mFaceMask.begin(), mFaceMask.end(), 0) == nfaces1);
  POLY_ASSERT(std::accumulate(mCellMask.begin(), mCellMask.end(), 0) == ncells1);

  cleanMesh(cellMap, faceMap, nodeMap);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
MeshEditor<Dimension, RealType>::
cleanMesh(std::vector<unsigned>& cellMap,
          std::vector<unsigned>& faceMap,
          std::vector<unsigned>& nodeMap) {

  // Pre-conditions
  POLY_ASSERT(!mCellMask.empty() and !mFaceMask.empty() and !mNodeMask.empty());

  // Original mesh sizes
  const unsigned ncells0 = mMesh.cells.size();
  const unsigned nfaces0 = mMesh.faces.size();
  const unsigned nnodes0 = mMesh.nodes.size()/Dimension;
  const unsigned ncells1 = std::accumulate(mCellMask.begin(), mCellMask.end(), 0);
  const unsigned nfaces1 = std::accumulate(mFaceMask.begin(), mFaceMask.end(), 0);
  const unsigned nnodes1 = std::accumulate(mNodeMask.begin(), mNodeMask.end(), 0);
  
  // // Count repeated node IDs to figure out which node positions need averaging
  // internal::CounterMap<unsigned> nodeCount;
  // for (unsigned i = 0; i != nnodes0; ++i)  ++nodeCount[nodeMap[i]];

  // // Reconstruct the nodes.
  // std::vector<RealType> newNodes(Dimension*nnodes1, 0.0);
  // for (unsigned i = 0; i != nnodes0; ++i) {
  //   const unsigned newID = nodeMap[i];
  //   POLY_ASSERT(nodeCount[newID] == 1 or nodeCount[newID] == 2);
  //   if (1) {//nodeCount[newID] == 1) {
  //     for (unsigned j = 0; j != Dimension; ++j) 
  //       newNodes[Dimension*newID+j] = mMesh.nodes[Dimension*i+j];
  //   } else {
  //     for (unsigned j = 0; j != Dimension; ++j)
  //       newNodes[Dimension*newID+j] += 0.5*mMesh.nodes[Dimension*i+j];
  //   }
  // }
  // mMesh.nodes = newNodes;
  
  std::vector<RealType> newNodes;
  newNodes.reserve(Dimension*nnodes1);
  for (unsigned i = 0; i != nnodes0; ++i) {
    if (mNodeMask[i] == 1) {
      std::copy(&mMesh.nodes[Dimension*i], &mMesh.nodes[Dimension*(i+1)], back_inserter(newNodes));
    }
  }
  mMesh.nodes = newNodes;

  // Reconstruct the faces.
  std::vector<std::vector<unsigned> > newFaces;
  std::vector<std::vector<int> > newFaceCells;
  newFaces.reserve(nfaces1);
  newFaceCells.reserve(nfaces1);
  for (unsigned i = 0; i != nfaces0; ++i) {
    if (mFaceMask[i] == 1) {
      newFaces.push_back(mMesh.faces[i]);
      for (std::vector<unsigned>::iterator itr = newFaces.back().begin();
           itr != newFaces.back().end(); 
           ++itr) {
         //POLY_ASSERT(nodeMap.find(*itr) != nodeMap.end());
        *itr = nodeMap[*itr];
      }
      POLY_ASSERT(mMesh.faceCells[i].size() == 1 or
		  mMesh.faceCells[i].size() == 2);
      newFaceCells.push_back(std::vector<int>());
      unsigned fc = (mMesh.faceCells[i][0] < 0 ? ~mMesh.faceCells[i][0] : mMesh.faceCells[i][0]);
      if (mCellMask[fc] == 1) newFaceCells.back().push_back(mMesh.faceCells[i][0] < 0 ?
                                                           ~cellMap[fc] :
                                                           cellMap[fc]);
      if (mMesh.faceCells[i].size() == 2) {
        fc = (mMesh.faceCells[i][1] < 0 ? ~mMesh.faceCells[i][1] : mMesh.faceCells[i][1]);
        if (mCellMask[fc] == 1) newFaceCells.back().push_back(mMesh.faceCells[i][1] < 0 ?
                                                              ~cellMap[fc] :
                                                              cellMap[fc]);
      }
      POLY_ASSERT(newFaceCells.back().size() == 1 or
		  newFaceCells.back().size() == 2);
    }
  }
  mMesh.faces = newFaces;
  mMesh.faceCells = newFaceCells;
  

  // Reconstruct the cells.  
  std::vector<std::vector<int> > newCells;
  newCells.resize(ncells1);
  for (unsigned i = 0; i != ncells0; ++i) {
    if (mCellMask[i] == 1) {
      for (std::vector<int>::iterator itr = mMesh.cells[i].begin();
           itr != mMesh.cells[i].end();
           ++itr) {
        const int iface = (*itr < 0 ? ~(*itr) : *itr);
        //POLY_ASSERT(faceMap.find(iface) != faceMap.end());
        if (mFaceMask[iface] == 1) {
          // int newFace = (*itr >= 0 ? faceMap[iface] : ~faceMap[iface]);
          newCells[i].push_back( *itr < 0 ? ~faceMap[iface] : faceMap[iface]);
        }
      }
    }
  }
  mMesh.cells = newCells;  

  // Update the shared nodes and faces.
  const unsigned numNeighbors = mMesh.sharedNodes.size();
  POLY_ASSERT(mMesh.sharedFaces.size() == numNeighbors);
  for (unsigned idomain = 0; idomain != numNeighbors; ++idomain) {
    std::vector<unsigned> newNodes, newFaces;
    for (std::vector<unsigned>::iterator itr = mMesh.sharedNodes[idomain].begin();
         itr != mMesh.sharedNodes[idomain].end();
         ++itr) {
      if (mNodeMask[*itr] == 1) newNodes.push_back(nodeMap[*itr]);
    }
    for (std::vector<unsigned>::iterator itr = mMesh.sharedFaces[idomain].begin();
         itr != mMesh.sharedFaces[idomain].end();
         ++itr) {
      if (mFaceMask[*itr] == 1) newFaces.push_back(faceMap[*itr]);
    }
    mMesh.sharedNodes[idomain] = newNodes;
    mMesh.sharedFaces[idomain] = newFaces;
  }

  // If there was a convex hull in the mesh, it's probably no longer valid.
  mMesh.convexHull = PLC<Dimension, RealType>();

  // Clear out the masks. We're done.
  mCellMask.clear();
  mFaceMask.clear();
  mNodeMask.clear();

  // Post-conditions.
  POLY_ASSERT(mMesh.nodes.size()     == Dimension*nnodes1);
  POLY_ASSERT(mMesh.faces.size()     == nfaces1          );
  POLY_ASSERT(mMesh.faceCells.size() == nfaces1          );
  POLY_ASSERT(mMesh.cells.size()     == ncells1          );
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
MeshEditor<Dimension, RealType>::
computeMasks() {
  // Pre-conditions
  POLY_ASSERT(!mCellMask.empty() or !mFaceMask.empty() or !mNodeMask.empty());
  
  // Have cells, want faces and nodes
  if (!mCellMask.empty()) {
    computeFaceAndNodeMasks();
  }
  // Have faces, want cells and nodes
  else if(!mFaceMask.empty()) {
    computeCellAndNodeMasks();
  }
  // Have nodes, want cells and faces
  else {
    POLY_ASSERT(!mNodeMask.empty());
    computeCellAndFaceMasks();
  }

  // Post-conditions
  POLY_ASSERT(!mCellMask.empty() and !mFaceMask.empty() and !mNodeMask.empty());
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
MeshEditor<Dimension, RealType>::
computeFaceAndNodeMasks() {
  // Pre-conditions.
  const unsigned ncells0 = mMesh.cells.size();
  const unsigned nfaces0 = mMesh.faces.size();
  const unsigned nnodes0 = mMesh.nodes.size()/Dimension;
  POLY_ASSERT(mCellMask.size() == ncells0);
  POLY_ASSERT(mFaceMask.empty() and mNodeMask.empty());
  POLY_ASSERT(ncells0 == 0 or *max_element(mCellMask.begin(), mCellMask.end()) == 1);

  // Create masks for the nodes and faces.
  mFaceMask = std::vector<unsigned>(nfaces0, 0);
  mNodeMask = std::vector<unsigned>(nnodes0, 0);
  for (unsigned icell = 0; icell != ncells0; ++icell) {
    if (mCellMask[icell] == 1) {
      for (std::vector<int>::const_iterator faceItr = mMesh.cells[icell].begin();
           faceItr != mMesh.cells[icell].end();
           ++faceItr) {
        const unsigned iface = (*faceItr >= 0 ? *faceItr : ~(*faceItr));
        POLY_ASSERT(iface < nfaces0);
        mFaceMask[iface] = 1;
        for (std::vector<unsigned>::const_iterator nodeItr = mMesh.faces[iface].begin();
             nodeItr != mMesh.faces[iface].end();
             ++nodeItr) {
          const unsigned inode = *nodeItr;
          POLY_ASSERT(inode < nnodes0);
          mNodeMask[inode] = 1;
        }
      }
    }
  }
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
MeshEditor<Dimension, RealType>::
computeCellAndNodeMasks() {
  // Pre-conditions.
  POLY_ASSERT(mCellMask.empty() and mNodeMask.empty());
  POLY_ASSERT2(false, "Not currently implemented");
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
MeshEditor<Dimension, RealType>::
computeCellAndFaceMasks() {
  // Pre-conditions.
  POLY_ASSERT(mCellMask.empty() and mFaceMask.empty());
  POLY_ASSERT2(false, "Not currently implemented");
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Static initializations
//------------------------------------------------------------------------------
template<> unsigned MeshEditor<2, double>::minEdgesPerFace = 1;
template<> unsigned MeshEditor<3, double>::minEdgesPerFace = 3;

//------------------------------------------------------------------------------
// Explicit instantiation
//------------------------------------------------------------------------------
template class MeshEditor<2, double>;
template class MeshEditor<3, double>;
  
}// end polytope namespace



/*
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
MeshEditor<Dimension, RealType>::
deleteBoundaryNodes(const std::vector<unsigned>& nodesToDelete) {
  POLY_ASSERT2(false, "This routine is currently not implemented");

  // Pre-conditions
  POLY_ASSERT(!nodesToDelete.empty());

  const unsigned nnodes0 = mMesh.nodes.size()/Dimension;
  const unsigned nfaces0 = mMesh.faces.size();
  const unsigned ncells0 = mMesh.cells.size();

  // Collect cell, faces, and nodes
  set<unsigned> boundaryCells, boundaryFaces, boundaryNodes;
  map<unsigned, set<unsigned> > boundaryNodesToBoundaryFaces;
  for (unsigned iface = 0; iface != mMesh.faces.size(); ++iface) {
    POLY_ASSERT(iface < mMesh.faceCells.size());
    POLY_ASSERT(mMesh.faceCells[iface].size() == 1 or
		mMesh.faceCells[iface].size() == 2);
    if (mMesh.faceCells[iface].size() == 1) {
      const unsigned icell = (mMesh.faceCells[iface][0] < 0) ?
	~mMesh.faceCells[iface][0] : mMesh.faceCells[iface][0];
      POLY_ASSERT(icell >= 0);
      POLY_ASSERT(icell < mMesh.cells.size());
      boundaryCells.insert(icell);
      boundaryFaces.insert(iface);
      POLY_ASSERT(!mMesh.faces[iface].empty());
      for (vector<unsigned>::const_iterator nodeItr = mMesh.faces[iface].begin(); 
	   nodeItr != mMesh.faces[iface].end(); ++nodeItr) {
	POLY_ASSERT(*nodeItr < nnodes0);
	boundaryNodes.insert(*nodeItr);
	boundaryNodesToBoundaryFaces[*nodeItr].insert(iface);
      }
    }
  }

  // Fill in the node mask
  mNodeMask = vector<unsigned>(nnodes0, 1);
  for (vector<unsigned>::const_iterator nodeItr = nodesToDelete.begin();
       nodeItr != nodesToDelete.end(); ++nodeItr) {
    POLY_ASSERT(*nodeItr < nnodes0);
    POLY_ASSERT2(boundaryNodes.find(*nodeItr) != boundaryNodes.end(),
		 "Node " << *nodeItr << " is flagged for deletion, but "
		 "it is NOT a boundary node!");
    mNodeMask[*nodeItr] = 0;
  }

  // Deleting nodes shouldn't remove any cells
  mCellMask = vector<unsigned>(ncells0, 1);

  // Deleting edges and faces is tricky. The 2D logic follows.
  // In 2D, edges=faces and a boundary node has at least 2 faces
  // These two faces correspond to boundary faces, each of which connects to a single boundary node
  // Collapse the node to the nearer of its two boundary node neighbors provided neither neighbor
  //    has already been flagged for deletion
  //
  // In 3D, edges!=faces.
  // A boundary node may have an arbitrary number of edges and faces, depending
  //    on the boundary
  // First gather all boundary node neighbors connecting to this node through
  //    boundary edges. (We obtain this list by looping over the nodes of its
  //    associated boundary faces and pick out the edges (i.e. node-node pairs)
  //    that also contain the node we want to delete.)
  // Collapse the node to the nearest of its neighboring boundary nodes,
  //    provided none of those nodes have also been flagged for deletion

}
//------------------------------------------------------------------------------
*/




/*
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
MeshEditor<Dimension, RealType>::
cleanBoundaryEdges(const RealType edgeTol) {
  // Pre-conditions
  POLY_ASSERT2(edgeTol > 0.0, "Specify a positive (non-zero) edge tolerance!");

  // flag edges smaller than the given tolerance in terms of face and node masks
  std::vector<unsigned> cellMap, faceMap, nodeMap;
  bool edgesClean = flagBoundaryEdgesForCleaning(edgeTol, cellMap, faceMap, nodeMap);
  if (!edgesClean) cleanMesh(cellMap, faceMap, nodeMap);

  // Post-conditions
  POLY_ASSERT(mCellMask.empty() and mFaceMask.empty() and mNodeMask.empty());
}
//------------------------------------------------------------------------------
*/



/*
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
bool
MeshEditor<Dimension, RealType>::
flagBoundaryEdgesForCleaning(const RealType edgeTol,
			     vector<unsigned>& cellMap,
			     vector<unsigned>& faceMap,
			     vector<unsigned>& nodeMap) {
  
  // Pre-conditions.
  POLY_ASSERT2(edgeTol > 0.0, "Specify a positive (non-zero) edge tolerance!");

  bool edgesClean = true;
  const unsigned ncells0 = mMesh.cells.size();
  const unsigned nfaces0 = mMesh.faces.size();
  const unsigned nnodes0 = mMesh.nodes.size()/Dimension;

  // Collect boundary cells
  set<unsigned> boundaryCells, boundaryNodes;
  for (unsigned iface = 0; iface != mMesh.faces.size(); ++iface) {
    POLY_ASSERT(iface < mMesh.faceCells.size());
    POLY_ASSERT(mMesh.faceCells[iface].size() == 1 or
		mMesh.faceCells[iface].size() == 2);
    if (mMesh.faceCells[iface].size() == 1) {
      const unsigned icell = (mMesh.faceCells[iface][0] < 0) ?
	~mMesh.faceCells[iface][0] : mMesh.faceCells[iface][0];
      boundaryCells.insert(icell);
      boundaryNodes.insert(mMesh.faces[iface].begin(),
			   mMesh.faces[iface].end());
    }
  }

  // Map edges to a mesh ID, compute maximum edge lengths per cell, 
  // and map cells to edges
  EdgeMap edgeToID;
  map<unsigned, vector<unsigned> > boundaryCellsToBoundaryEdges;
  vector<RealType> edgeLength;
  RealType length;
  unsigned inode0, inode1;
  int edgeCount = 0;
  for (set<unsigned>::const_iterator cellItr = boundaryCells.begin();
       cellItr != boundaryCells.end(); ++cellItr) {
    RealType maxEdgeLength = 0.0;
    vector<unsigned> cellEdgeIDs;
    for (vector<int>::const_iterator itr = mMesh.cells[*cellItr].begin();
	 itr != mMesh.cells[*cellItr].end(); ++itr) {
      const unsigned iface = (*itr < 0) ? ~(*itr) : *itr;
      POLY_ASSERT(iface < mMesh.faces.size());
      const unsigned nfaceNodes = mMesh.faces[iface].size();
      POLY_ASSERT(nfaceNodes >= 2);
      const unsigned maxNodeIndex = (nfaceNodes == 2) ? 1 : nfaceNodes;
      for (unsigned inode = 0; inode != maxNodeIndex; ++inode) {
        inode0 = mMesh.faces[iface][inode];
        inode1 = mMesh.faces[iface][(inode+1) % nfaceNodes];
        POLY_ASSERT(inode0 != inode1);
	POLY_ASSERT(inode0 < mMesh.nodes.size()/Dimension and
                    inode1 < mMesh.nodes.size()/Dimension);
        const EdgeType edge = internal::hashEdge(inode0, inode1);
	typename EdgeMap::left_const_iterator lItr = edgeToID.left.find(edge);
	if (lItr == edgeToID.left.end()) {
	  edgeToID.insert(typename EdgeMap::value_type(edge, edgeCount));
	  length = geometry::distance<Dimension, RealType>
	    (&mMesh.nodes[Dimension*inode0], &mMesh.nodes[Dimension*inode1]);
	  POLY_ASSERT(length > 0.0);
	  maxEdgeLength = std::max(maxEdgeLength, length);
	  POLY_ASSERT(iface < mMesh.faceCells.size());
	  if (mMesh.faceCells[iface].size() == 1) {
	    boundaryCellsToBoundaryEdges[*cellItr].push_back(edgeCount);
	  }
	  //edgeLength.push_back(length);
	  //cellEdgeIDs.push_back(edgeCount);
	  ++edgeCount;
	} else {
	  cellEdgeIDs.push_back(lItr->second);
	}
      }
    }
    cellToEdges[*cellItr] = cellEdgeIDs;
  }
  POLY_ASSERT(edgeCount == edgeLength.size());
  POLY_ASSERT(cellToEdges.size() == ncells0);
  const unsigned nedges0 = edgeCount;



  //
        typename EdgeMap::left_const_iterator lItr = edgeToID.left.find(edge);
	if (lItr == edgeToID.left.end()) {
	  edgeToID.insert(typename EdgeMap::value_type(edge, edgeCount));
	  length = geometry::distance<Dimension, RealType>
	    (&mMesh.nodes[Dimension*inode0], &mMesh.nodes[Dimension*inode1]);
	  edgeLength.push_back(length);
	  cellEdgeIDs.push_back(edgeCount);
	  ++edgeCount;
	} else {
	  cellEdgeIDs.push_back(lItr->second);
	}
      }
    }
    cellToEdges[icell] = cellEdgeIDs;
  }
  POLY_ASSERT(edgeCount == edgeLength.size());
  POLY_ASSERT(cellToEdges.size() == ncells0);
  const unsigned nedges0 = edgeCount;

  //





  // Compute the maximum edge length for the cells around an edge
  vector<RealType> maxCellEdgeLength(nedges0, 0.0);
  for (set<unsigned>::const_iterator cellItr = boundaryCells.begin();
       cellItr != boundaryCells.end(); ++cellItr) {
    RealType maxCellEdge = 0.0;
    vector<unsigned> edgeIDs = cellToEdges[*cellItr];
    for (vector<unsigned>::const_iterator itr = edgeIDs.begin();
	 itr != edgeIDs.end();
	 ++itr) maxCellEdge = max(maxCellEdge, edgeLength[*itr]);
    for (vector<unsigned>::const_iterator itr = edgeIDs.begin();
	 itr != edgeIDs.end();
	 ++itr) maxCellEdgeLength[*itr] = max(maxCellEdgeLength[*itr], maxCellEdge);
  }
  
  // Flag the edges we want to remove.
  vector<unsigned> edgeMask(nedges0, 1);
  vector<unsigned> nodeCollapse(nnodes0);
  EdgeType edge;
  mNodeMask = vector<unsigned>(nnodes0, 1);
  for (unsigned i = 0; i != nnodes0; ++i) nodeCollapse[i] = i;
  for (unsigned iedge = 0; iedge != nedges0; ++iedge) {
    typename EdgeMap::right_const_iterator itr = edgeToID.right.find(iedge);
    POLY_ASSERT(itr != edgeToID.right.end());
    edge = itr->second;
    inode0 = edge.first;
    inode1 = edge.second;
    if (edgeLength[iedge] < edgeTol*maxCellEdgeLength[iedge] and
	mNodeMask[inode0] == 1 and mNodeMask[inode1] == 1) {
      edgesClean = false;
      edgeMask[iedge] = 0;
      mNodeMask[inode0] = 2;
      mNodeMask[inode1] = 0;
      nodeCollapse[inode1] = inode0;
    }
  }
  replace_if(mNodeMask.begin(), mNodeMask.end(), bind2nd(equal_to<unsigned>(), 2), 1);

  
#ifdef POLYTOPE_ENABLE_MPI
  // Parallel configuration.
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // TODO: Put in an all-reduced switch based on the edgesClean flag. After all, if
  //       all processors have clean mesh, why bother communicating?

  // Construct a list of shared edges, sorted the same on each proc. Construct and 
  // communicate the masks and make the local edge mask consistent across domains.
  // Of course, none of this is necessary if we're running in serial, so...
  if (numProcs > 1) {

    // // Create reverse loopup from edge to collection of neighboring cells
    // map<EdgeType, vector<unsigned> > edgeToCells;
    // for (map<unsigned, vector<unsigned> >::const_iterator cellItr = cellToEdges.begin();
    //      cellItr != cellToEdges.end(); ++cellItr) {
    //   vector<unsigned> 
    // }

    // Construct shared edge list
    vector<vector<EdgeType> > sharedEdges   (mMesh.neighborDomains.size());
    vector<vector<unsigned> > sharedEdgeMask(mMesh.neighborDomains.size());
    vector<vector<unsigned> > sharedEdgeReverse(nedges0);
    vector<vector<unsigned> > sharedNodeReverse(nnodes0);
    vector<map<unsigned,unsigned> > domainNodeMaps(mMesh.neighborDomains.size());
    unsigned ineighbor = 0;
    for (vector<unsigned>::const_iterator otherItr = mMesh.neighborDomains.begin();
         otherItr != mMesh.neighborDomains.end(); ++otherItr, ++ineighbor) {
      
      // Create reverse lookup from nodeID to position in sharedNodes list
      map<unsigned, unsigned> nodeIndexToSharePosition;
      for (vector<unsigned>::const_iterator nodeItr = mMesh.sharedNodes[ineighbor].begin();
           nodeItr != mMesh.sharedNodes[ineighbor].end(); ++nodeItr) {
        POLY_ASSERT(nodeIndexToSharePosition.find(*nodeItr) == nodeIndexToSharePosition.end());
        unsigned result = nodeIndexToSharePosition.size();
        nodeIndexToSharePosition[*nodeItr] = result;
        sharedNodeReverse[*nodeItr].push_back(ineighbor);
      }
      POLY_ASSERT(nodeIndexToSharePosition.size() == mMesh.sharedNodes[ineighbor].size());
      domainNodeMaps[ineighbor] = nodeIndexToSharePosition;

      // Consider all edges on this domain. The shared ones consist of two shared nodes.
      for (typename EdgeMap::left_const_iterator edgeItr = edgeToID.left.begin();
           edgeItr != edgeToID.left.end(); ++edgeItr) {
        const EdgeType edge  = edgeItr->first;
        const unsigned iedge = edgeItr->second;
        inode0 = edge.first;
        inode1 = edge.second;
        POLY_ASSERT(inode0 < inode1);
        POLY_ASSERT(inode0 < mMesh.nodes.size()/Dimension and
                    inode1 < mMesh.nodes.size()/Dimension);
        map<unsigned,unsigned>::const_iterator it0 = nodeIndexToSharePosition.find(inode0);
        map<unsigned,unsigned>::const_iterator it1 = nodeIndexToSharePosition.find(inode1);
        if (it0 != nodeIndexToSharePosition.end() and 
            it1 != nodeIndexToSharePosition.end()) {
          const EdgeType sharedEdge = internal::hashEdge(it0->second, it1->second);
          sharedEdges[ineighbor].push_back(sharedEdge);
          sharedEdgeReverse[iedge].push_back(ineighbor);
        }
      }

      // Sort the sharedEdges list
      sort(sharedEdges[ineighbor].begin(), sharedEdges[ineighbor].end(),
           internal::pairCompare<unsigned, unsigned>);        

      // Now make the corresponding mask for each sharedEdge list
      sharedEdgeMask[ineighbor].resize(sharedEdges[ineighbor].size());
      unsigned edgeCount = 0;
      for (vector<EdgeType>::const_iterator sEdgeItr = sharedEdges[ineighbor].begin();
           sEdgeItr != sharedEdges[ineighbor].end(); ++sEdgeItr, ++edgeCount) {
        const unsigned ind0 = sEdgeItr->first;
        const unsigned ind1 = sEdgeItr->second;
        POLY_ASSERT(ind0 < mMesh.sharedNodes[ineighbor].size() and
                    ind1 < mMesh.sharedNodes[ineighbor].size());
        inode0 = mMesh.sharedNodes[ineighbor][ind0];
        inode1 = mMesh.sharedNodes[ineighbor][ind1];
        POLY_ASSERT(inode0 < mMesh.nodes.size()/Dimension and
                    inode1 < mMesh.nodes.size()/Dimension);
        const EdgeType edge = internal::hashEdge(inode0, inode1);
        typename EdgeMap::left_const_iterator itr = edgeToID.left.find(edge);
        POLY_ASSERT(itr != edgeToID.left.end());
        const unsigned iedge = itr->second;
        POLY_ASSERT(iedge < nedges0);
        sharedEdgeMask[ineighbor][edgeCount] = edgeMask[iedge];
      }
      POLY_ASSERT(edgeCount == sharedEdges[ineighbor].size());
    }
    
    // Send the sharedEdgeMask to each neighbor asynchronously
    vector<char> localBuffer;
    vector<MPI_Request> sendRequests;
    sendRequests.reserve(2*mMesh.neighborDomains.size());
    ineighbor = 0;
    for (vector<unsigned>::const_iterator otherItr = mMesh.neighborDomains.begin();
         otherItr != mMesh.neighborDomains.end(); ++otherItr, ++ineighbor) {
      serialize(sharedEdgeMask[ineighbor], localBuffer);
      unsigned localBufferSize = localBuffer.size();
      const unsigned otherProc = *otherItr;
      sendRequests.push_back(MPI_Request());
      MPI_Isend(&localBufferSize, 1, MPI_UNSIGNED, otherProc, 1, MPI_COMM_WORLD, &sendRequests.back());
      if (localBufferSize > 0) {
        sendRequests.push_back(MPI_Request());
        MPI_Isend(&localBuffer.front(), localBufferSize, MPI_CHAR, otherProc, 2, MPI_COMM_WORLD, &sendRequests.back());
      }
    }
    POLY_ASSERT(sendRequests.size() <= 2*mMesh.neighborDomains.size());

    // Get the masks from the neighbors
    vector<vector<unsigned> > neighborEdgeMasks(mMesh.neighborDomains.size());
    ineighbor = 0;
    for (vector<unsigned>::const_iterator otherItr = mMesh.neighborDomains.begin();
         otherItr != mMesh.neighborDomains.end(); ++otherItr, ++ineighbor) {
      const unsigned otherProc = *otherItr;
      unsigned bufSize;
      MPI_Status recvStatus1, recvStatus2;
      MPI_Recv(&bufSize, 1, MPI_UNSIGNED, otherProc, 1, MPI_COMM_WORLD, &recvStatus1);
      if (bufSize > 0) {
        vector<char> buffer(bufSize, '\0');
        MPI_Recv(&buffer.front(), bufSize, MPI_CHAR, otherProc, 2, MPI_COMM_WORLD, &recvStatus2);
        vector<char>::const_iterator itr = buffer.begin();
        vector<unsigned> otherEdgeMask;
        deserialize(otherEdgeMask, itr, buffer.end());
        POLY_ASSERT(itr == buffer.end());
        neighborEdgeMasks[ineighbor] = otherEdgeMask;
      }
    }
    
    // Make sure all the sends are complete
    vector<MPI_Status> sendStatus(sendRequests.size());
    MPI_Waitall(sendRequests.size(), &sendRequests.front(), &sendStatus.front());
    
    // Make the edge masks consistent across domain boundaries. If one domain 
    // says delete (mask=0) and one says keep (mask=1), we always keep.
    for (unsigned i = 0; i != mMesh.neighborDomains.size(); ++i) {
      unsigned edgeCount = 0;
      for (vector<EdgeType>::const_iterator sEdgeItr = sharedEdges[i].begin();
           sEdgeItr != sharedEdges[i].end(); ++sEdgeItr, ++edgeCount) {
        const unsigned ind0 = sEdgeItr->first;
        const unsigned ind1 = sEdgeItr->second;
        POLY_ASSERT(ind0 < mMesh.sharedNodes[i].size() and
                    ind1 < mMesh.sharedNodes[i].size());
        inode0 = mMesh.sharedNodes[i][ind0];
        inode1 = mMesh.sharedNodes[i][ind1];
        POLY_ASSERT(inode0 < mMesh.nodes.size()/Dimension and
                    inode1 < mMesh.nodes.size()/Dimension);
        const EdgeType edge = internal::hashEdge(inode0, inode1);
        typename EdgeMap::left_const_iterator itr = edgeToID.left.find(edge);
        POLY_ASSERT(itr != edgeToID.left.end());
        const unsigned iedge = itr->second;
        POLY_ASSERT(iedge < nedges0);
        //POLY_ASSERT(edgeCount < neighborEdgeMasks[i].size());
        edgeMask[iedge] = max(edgeMask[iedge], neighborEdgeMasks[i][edgeCount]);
      }
      POLY_ASSERT(edgeCount == sharedEdges[i].size());
    }

    // Ammend the edgesClean flag if we're not cleaning after all.
    if (std::accumulate(edgeMask.begin(), edgeMask.end(), 0) == nedges0)
       edgesClean = true;

    // We've invalidated the nodeCollapse and nodeMask lists. Make them consistent.
    mNodeMask = vector<unsigned>(nnodes0, 1);
    if (!edgesClean) {
      for (unsigned i = 0; i != nnodes0; ++i) nodeCollapse[i] = i;
      for (unsigned iedge = 0; iedge != nedges0; ++iedge) {
        typename EdgeMap::right_const_iterator itr = edgeToID.right.find(iedge);
        POLY_ASSERT(itr != edgeToID.right.end());
        if (edgeMask[iedge] == 0) {
          const EdgeType edge = itr->second;
          const bool isShared = (sharedEdgeReverse[iedge].size() == 0) ? false : true;
          inode0 = edge.first;
          inode1 = edge.second;
          const bool node0Shared = (sharedNodeReverse[inode0].size() == 0) ? false : true;
          const bool node1Shared = (sharedNodeReverse[inode1].size() == 0) ? false : true;

          if (isShared) {
            POLY_ASSERT2(sharedEdgeReverse[iedge].size() == 1, "I think this logic only"
                         << " holds for 2D meshes. An edge is a face and can be shared by"
                         << " at most two domains. In 3D, 3 or more domains can share an"
                         << " edge. It's not clear if the orderings of the shared node"
                         << " lists are consistent among 3 or more domains");
            const unsigned ineighbor = sharedEdgeReverse[iedge][0];
            const unsigned ind0 = domainNodeMaps[ineighbor][inode0];
            const unsigned ind1 = domainNodeMaps[ineighbor][inode1];
            if (ind1 < ind0) {inode0 = edge.second; inode1 = edge.first;}
          }

          else if (node1Shared and !node0Shared) {
             inode0 = edge.second; inode1 = edge.first;
          }

          mNodeMask[inode1] = 0;
          nodeCollapse[inode1] = inode0;
        }
      }
    }

  } // Done with the multi-processor stuff
#endif

        
  // Check if there are any faces that need to be removed (for 3D meshes).
  mFaceMask = std::vector<unsigned>(nfaces0, 1);
  if (!edgesClean) {
    unsigned numActiveEdges;
    for (unsigned iface = 0; iface != nfaces0; ++iface) {
      numActiveEdges = 0;
      const unsigned maxNodeIndex = (mMesh.faces[iface].size() == 2) ? 1 : mMesh.faces[iface].size();
      for (unsigned inode = 0; inode != maxNodeIndex; ++inode) {
        inode0 = mMesh.faces[iface][inode];
        inode1 = mMesh.faces[iface][(inode+1) % mMesh.faces[iface].size()];
        const EdgeType edge = internal::hashEdge(inode0, inode1);
	typename EdgeMap::left_const_iterator itr = edgeToID.left.find(edge);
	POLY_ASSERT(itr != edgeToID.left.end());
	numActiveEdges += edgeMask[itr->second];
      }
      if( numActiveEdges < minEdgesPerFace) mFaceMask[iface] = 0;
    }
  }
  
  // Build up the node map from original indices to new indices
  unsigned nodeCount = 0;
  nodeMap.resize(nnodes0);
  for (unsigned i = 0; i != nnodes0; ++i) {
    nodeMap[i] = nodeCount;
    nodeCount += mNodeMask[i];
  }
  for (unsigned i = 0; i != nnodes0; ++i) {
    nodeMap[i] = nodeMap[nodeCollapse[i]];
  }
  POLY_ASSERT(nodeCount == std::accumulate(mNodeMask.begin(), mNodeMask.end(), 0));

  // Build up the face map from original indices to new indices
  faceMap.resize(nfaces0);
  unsigned faceCount = 0;
  for (unsigned iface = 0; iface != nfaces0; ++iface) {
    faceMap[iface] = faceCount;
    faceCount += mFaceMask[iface];
  }
  POLY_ASSERT(faceCount <= nfaces0);

  // Deleting edges (hopefully) does not result in a deleted cell
  mCellMask = std::vector<unsigned>(ncells0, 1);
  cellMap.resize(ncells0);
  for (unsigned i = 0; i != ncells0; ++i) cellMap[i] = i;

  //Post-conditions
  if (edgesClean) {
     POLY_ASSERT(std::accumulate(mNodeMask.begin(), mNodeMask.end(), 0) == nnodes0);
     POLY_ASSERT(std::accumulate(mFaceMask.begin(), mFaceMask.end(), 0) == nfaces0);
     POLY_ASSERT(std::accumulate(mCellMask.begin(), mCellMask.end(), 0) == ncells0);
     mNodeMask.clear();  mFaceMask.clear();  mCellMask.clear();
  } else {
    POLY_ASSERT(std::accumulate(mNodeMask.begin(), mNodeMask.end(), 0) < nnodes0);
  }
  
  return edgesClean;
}
//------------------------------------------------------------------------------
*/
