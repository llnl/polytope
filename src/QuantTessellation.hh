//------------------------------------------------------------------------------
// An internal handy intermediate representation of a tessellation.
//------------------------------------------------------------------------------
#ifndef __Polytope_QuantTessellation__
#define __Polytope_QuantTessellation__

#include <algorithm>
#include "polytope_internal.hh"
#include "DimensionTraits.hh"

namespace polytope {
namespace internal {

template<int Dimension, typename RealType>
class QuantTessellation {
public:
  typedef typename DimensionTraits<Dimension, RealType>::CoordHash CoordHash;
  typedef typename DimensionTraits<Dimension, RealType>::CoordHash PointHash;
  typedef typename DimensionTraits<Dimension, RealType>::IntPoint  IntPoint;
  typedef typename DimensionTraits<Dimension, RealType>::RealPoint RealPoint;
  typedef std::pair<int, int> EdgeHash;
  typedef std::vector<unsigned> FaceHash;

  // The normalized generator coordinates.
  std::vector<RealType> generators;

  // The bounds for hashing positions.
  RealPoint low_labframe, high_labframe;
  RealPoint low_inner, high_inner, low_outer, high_outer;

  // The degeneracy we're using for quantizing.
  RealType degeneracy;

  //----------------------------------------------------------------------------
  // The mesh elements and connectivity.
  //----------------------------------------------------------------------------
  std::map<PointHash, int> point2id;              // PointHash -> unique point ID
  std::map<EdgeHash, int> edge2id;                // EdgeHash  -> unique edge ID
  std::map<FaceHash, int> face2id;                // FaceHash  -> unique face ID
  std::vector<PointHash> points;                  // Hashed node positions.
  std::vector<EdgeHash> edges;                    // Hashed edges (node index pairs).
  std::vector<std::vector<int> > faces;           // Faces made of edges (with orientation)
  std::vector<std::vector<int> > cells;           // Cells made of faces (with orientation)
  std::vector<unsigned> infNodes;                 // Indices of nodes projected to the infSphere
  std::vector<unsigned> infEdges;                 // Indices of edges projected to the infSphere
  std::vector<unsigned> infFaces;                 // Indices of faces projected to the infSphere

  //----------------------------------------------------------------------------
  // Hash the given position.
  //----------------------------------------------------------------------------
  PointHash hashPosition(const RealPoint& p) const {
    return geometry::Hasher<Dimension, RealType>::hashPosition(const_cast<RealType*>(&(p.x)), 
                                                               const_cast<RealType*>(&low_inner.x), const_cast<RealType*>(&high_inner.x), 
                                                               const_cast<RealType*>(&low_outer.x), const_cast<RealType*>(&high_outer.x),
                                                               degeneracy);
  }
  RealPoint unhashPosition(const PointHash ip) const {
    RealPoint result;
    geometry::Hasher<Dimension, RealType>::unhashPosition(&result.x, 
                                                          const_cast<RealType*>(&low_inner.x), const_cast<RealType*>(&high_inner.x), 
                                                          const_cast<RealType*>(&low_outer.x), const_cast<RealType*>(&high_outer.x), 
                                                          ip, degeneracy);
    return result;
  }
  IntPoint hashedPosition(const PointHash ip) const {
    IntPoint result;
    geometry::Hasher<Dimension, RealType>::hashedPosition(&result.x, ip);
    return result;
  }

  //----------------------------------------------------------------------------
  // Add new elements, and return the unique index.
  //----------------------------------------------------------------------------
  // Nodes
  int addNewNode(const PointHash ix) {
    const int k = point2id.size();
    const int result = internal::addKeyToMap(ix, point2id);
    if (result == k) {
      points.push_back(ix);
    }
    POLY_ASSERT(points.size() == point2id.size());
    return result;
  }

  int addNewNode(const RealPoint& x) {
    return this->addNewNode(hashPosition(x));
  }

  // Edges.
  int addNewEdge(const EdgeHash& x) {
    const int k = edge2id.size();
    const int result = internal::addKeyToMap(x, edge2id);
    if (result == k) {
      edges.push_back(x);
    }
    POLY_ASSERT(edges.size() == edge2id.size());
    return result;
  }

  // Faces.
  // Note this is a little different than above.  A FaceHash is not the same
  // as the signed, oriented, and ordered list of edges that constitute a face!
  // Rather a FaceHash is the sorted positive IDs of the edges in the face.
  int addNewFace(const std::vector<int>& x) {
    FaceHash fhash;
    for (unsigned i = 0; i != x.size(); ++i) fhash.push_back(x[i] < 0 ? ~x[i] : x[i]);
    std::sort(fhash.begin(), fhash.end());
    const int k = face2id.size();
    const int result = internal::addKeyToMap(fhash, face2id);
    if (result == k) {
      faces.push_back(x);
    }
    POLY_ASSERT(faces.size() == face2id.size());
    return result;
  }

  //----------------------------------------------------------------------------
  // Floating position for a point (normalized coordinates).
  //----------------------------------------------------------------------------
  RealPoint nodePosition(const unsigned i) {
    POLY_ASSERT(i < points.size());
    return unhashPosition(points[i]);
  }

//   //----------------------------------------------------------------------------
//   // Integer position for a point (normalized coordinates).
//   //----------------------------------------------------------------------------
//   IntPoint nodePosition(const unsigned i) {
//     POLY_ASSERT(i < points.size());
//     return hashedPosition(points[i]);
//   }

  //----------------------------------------------------------------------------
  // Floating position for a point (lab frame).
  //----------------------------------------------------------------------------
  RealPoint labNodePosition(const unsigned i) {
    POLY_ASSERT(i < points.size());
    RealPoint result = nodePosition(i);
    for (unsigned j = 0; j != Dimension; ++j) {
      result[j] = result[j]*(high_labframe[j] - low_labframe[j]) + low_labframe[j];
    }
    return result;
  }

  //----------------------------------------------------------------------------
  // Floating position for a point (lab frame).
  //----------------------------------------------------------------------------
  RealPoint labNodePositionCollinear(const unsigned i) {
    POLY_ASSERT(i < points.size());
    RealPoint result = nodePosition(i);
    for (unsigned j = 0; j != Dimension; ++j) {
      if (low_labframe[j] == high_labframe[j]) {
	result[j] = result[j] + low_labframe[j];
      } else {
	result[j] = result[j]*(high_labframe[j] - low_labframe[j]) + low_labframe[j];
      }
    }
    return result;
  }

  //----------------------------------------------------------------------------
  // Floating position for an edge.
  //----------------------------------------------------------------------------
  RealPoint edgePosition(const EdgeHash& ehash) {
    RealPoint result = nodePosition(ehash.first) + nodePosition(ehash.second);
    result *= 0.5;
    return result;
  }

  //----------------------------------------------------------------------------
  // Compute the node->edge connectivity.
  //----------------------------------------------------------------------------
  std::vector<std::vector<unsigned> > nodeEdges() const {
    std::vector<std::vector<unsigned> > result(points.size());
    for (unsigned iedge = 0; iedge != edges.size(); ++iedge) {
      POLY_ASSERT(edges[iedge].first < points.size());
      POLY_ASSERT(edges[iedge].second < points.size());
      result[edges[iedge].first].push_back(iedge);
      result[edges[iedge].second].push_back(iedge);
    }
    return result;
  }

  //----------------------------------------------------------------------------
  // Compute the edge->face connectivity.
  //----------------------------------------------------------------------------
  std::vector<std::vector<int> > edgeFaces() const {
    std::vector<std::vector<int> > result(edges.size());
    for (int iface = 0; iface != faces.size(); ++iface) {
      for (unsigned k = 0; k != faces[iface].size(); ++k) {
        if (faces[iface][k] < 0) {
          POLY_ASSERT(~faces[iface][k] < edges.size());
          result[~faces[iface][k]].push_back(~iface);
        } else {
          POLY_ASSERT(faces[iface][k] < edges.size());
          result[faces[iface][k]].push_back(iface);
        }
      }
    }
    return result;
  }

  //----------------------------------------------------------------------------
  // Compute the face->cell connectivity.
  //----------------------------------------------------------------------------
  std::vector<std::vector<int> > faceCells() const {
    std::vector<std::vector<int> > result(faces.size());
    for (int icell = 0; icell != cells.size(); ++icell) {
      for (unsigned k = 0; k != cells[icell].size(); ++k) {
        if (cells[icell][k] < 0) {
          POLY_ASSERT(~cells[icell][k] < faces.size());
          result[~cells[icell][k]].push_back(~icell);
        } else {
          POLY_ASSERT(cells[icell][k] < faces.size());
          result[cells[icell][k]].push_back(icell);
        }
      }
    }
    return result;
  }

  //----------------------------------------------------------------------------
  // Convert our internal data to a standard polytope Tessellation.
  //----------------------------------------------------------------------------
  void tessellation(Tessellation<Dimension, RealType>& mesh) {
    mesh.clear();

    // Nodes.
    mesh.nodes.resize(Dimension*points.size());
    if (low_labframe.x == high_labframe.x or
	low_labframe.y == high_labframe.y) {
      for (unsigned i = 0; i != points.size(); ++i) {
	RealPoint p = labNodePositionCollinear(i);
	std::copy(&p.x, &p.x + Dimension, &mesh.nodes[Dimension*i]);
      }
    } else {
      for (unsigned i = 0; i != points.size(); ++i) {
	RealPoint p = labNodePosition(i);
	std::copy(&p.x, &p.x + Dimension, &mesh.nodes[Dimension*i]);
      }
    }

    // Faces.
    mesh.faces.reserve(faces.size());
    for (unsigned i = 0; i != faces.size(); ++i) {
      mesh.faces.push_back(std::vector<unsigned>());
      if (faces[i].size() == 1) {
	mesh.faces[i].push_back(edges[faces[i][0]].first);
	mesh.faces[i].push_back(edges[faces[i][0]].second);
      } else {
	for (unsigned j = 0; j != faces[i].size(); ++j) {
	  if (faces[i][j] >= 0) {
	    mesh.faces[i].push_back(edges[faces[i][j]].first);
	  } else {
	    mesh.faces[i].push_back(edges[~faces[i][j]].second);
	  }
	  POLY_ASSERT(mesh.faces[i].back() < mesh.nodes.size()/Dimension);
	}
	POLY_ASSERT(mesh.faces[i].size() == faces[i].size());
      }
    }

    // Much of our data can simply be copied over wholesale.
    mesh.cells = cells;
    mesh.infNodes = infNodes;
    mesh.infFaces = infFaces;
    mesh.faceCells = this->faceCells();
  }

  //----------------------------------------------------------------------------
  // A contract heavy validity check.
  //----------------------------------------------------------------------------
  void assertValid() const {
    POLY_BEGIN_CONTRACT_SCOPE;
    {
      if (Dimension == 3) {
        const QuantTessellation<Dimension, RealType> qmesh = *this;
        const unsigned numGenerators = qmesh.generators.size()/Dimension;
        const std::vector<std::vector<unsigned> > nodeEdges = qmesh.nodeEdges();
        const std::vector<std::vector<int> > edgeFaces = qmesh.edgeFaces();
        const std::vector<std::vector<int> > faceCells = qmesh.faceCells();
        POLY_ASSERT(qmesh.points.size() == qmesh.point2id.size());
        POLY_ASSERT(qmesh.edges.size() == qmesh.edge2id.size());
        POLY_ASSERT(qmesh.faces.size() == qmesh.face2id.size());
        POLY_ASSERT(qmesh.cells.size() == numGenerators);
        POLY_ASSERT(nodeEdges.size() == qmesh.point2id.size());
        POLY_ASSERT(edgeFaces.size() == qmesh.edges.size());
        POLY_ASSERT(faceCells.size() == qmesh.faces.size());
        for (int i = 0; i != qmesh.points.size(); ++i) {
          for (int j = 0; j != nodeEdges[i].size(); ++j) {
            POLY_ASSERT(qmesh.edges[nodeEdges[i][j]].first == i or qmesh.edges[nodeEdges[i][j]].second == i);
          }
        }
        for (int i = 0; i != qmesh.edges.size(); ++i) {
          POLY_ASSERT(qmesh.edges[i].first < qmesh.points.size());
          POLY_ASSERT(qmesh.edges[i].second < qmesh.points.size());
          POLY_ASSERT(count(nodeEdges[qmesh.edges[i].first].begin(),
                            nodeEdges[qmesh.edges[i].first].end(), 
                            i) == 1);
          POLY_ASSERT(count(nodeEdges[qmesh.edges[i].second].begin(),
                            nodeEdges[qmesh.edges[i].second].end(), 
                            i) == 1);
          for (int j = 0; j != edgeFaces[i].size(); ++j) {
            const int iface = edgeFaces[i][j];
            POLY_ASSERT(internal::positiveID(iface) < qmesh.faces.size());
            if (iface < 0) {
              POLY_ASSERT(count(qmesh.faces[~iface].begin(), qmesh.faces[~iface].end(), ~i) == 1);
            } else {
              POLY_ASSERT(count(qmesh.faces[iface].begin(), qmesh.faces[iface].end(), i) == 1);
            }
          }
        }
        for (int i = 0; i != qmesh.faces.size(); ++i) {
          const unsigned nedges = qmesh.faces[i].size();
          POLY_ASSERT(nedges >= 3);
          for (int j = 0; j != nedges; ++j) {
            const int k = (j + 1) % nedges;
            const int iedge1 = qmesh.faces[i][j];
            const int iedge2 = qmesh.faces[i][k];
            POLY_ASSERT(internal::positiveID(iedge1) < qmesh.edges.size());
            POLY_ASSERT(internal::positiveID(iedge2) < qmesh.edges.size());
            if (iedge1 >= 0 and iedge2 >= 0) {
              POLY_ASSERT(qmesh.edges[iedge1].second == qmesh.edges[iedge2].first);
            } else if (iedge1 >= 0 and iedge2 < 0) {
              POLY_ASSERT(qmesh.edges[iedge1].second == qmesh.edges[~iedge2].second);
            } else if (iedge1 < 0 and iedge2 >= 0) {
              POLY_ASSERT(qmesh.edges[~iedge1].first == qmesh.edges[iedge2].first);
            } else {
              POLY_ASSERT(qmesh.edges[~iedge1].first == qmesh.edges[~iedge2].second);
            }
            if (iedge1 < 0) {
              POLY_ASSERT(count(edgeFaces[~iedge1].begin(), edgeFaces[~iedge1].end(), ~i) == 1);
            } else {
              POLY_ASSERT(count(edgeFaces[iedge1].begin(), edgeFaces[iedge1].end(), i) == 1);
            }
          }
          POLY_ASSERT(faceCells[i].size() == 1 or faceCells[i].size() == 2);
          for (int j = 0; j != faceCells[i].size(); ++ j) {
            const int icell = faceCells[i][j];
            POLY_ASSERT(internal::positiveID(icell) < numGenerators);
            if (icell < 0) {
              POLY_ASSERT(count(qmesh.cells[~icell].begin(), qmesh.cells[~icell].begin(), ~i) == 1);
            } else {
              POLY_ASSERT(count(qmesh.cells[icell].begin(), qmesh.cells[icell].begin(), i) == 1);
            }
          }
        }
        for (int i = 0; i != numGenerators; ++i) {
          const unsigned nfaces = qmesh.cells[i].size();
          POLY_ASSERT(nfaces >= 4);
          for (int j = 0; j != nfaces; ++j) {
            const int iface = qmesh.cells[i][j];
            POLY_ASSERT(internal::positiveID(iface) < qmesh.faces.size());
            if (iface < 0) {
              POLY_ASSERT(count(faceCells[~iface].begin(), faceCells[~iface].end(), ~i) == 1);
            } else {
              POLY_ASSERT(count(faceCells[iface].begin(), faceCells[iface].end(), i) == 1);
            }
          }
        }
      }
    }
    POLY_END_CONTRACT_SCOPE;
  }
};

}
}

#endif
