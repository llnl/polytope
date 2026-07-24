#ifndef __Polytope_Tessellation__
#define __Polytope_Tessellation__

#include <vector>
#include <set>
#include <iostream>
#include "PLC.hh"
#include "polytope_internal.hh"
#include "Cell.hh"
namespace polytope {

//! \class Mesh - A basic descriptor class for a topologically-consistent
//! arbitrary poly(gonal/hedral) mesh.
template<int Dimension, typename RealType>
class Tessellation {
  public:
  using RealCell = typename Cell<Dimension, RealType>::CellType;
  // Default constructor.
  Tessellation():
    points(),
    nodes(),
    cells(),
    faces(),
    infNodes(),
    boundaryNodes(),
    infFaces(),
    boundaryFaces(),
    faceCells(),
    convexHull() {}

  // Destructor.
  virtual ~Tessellation() {};

  //! Clears the tessellation, emptying it of all data.
  virtual void clear()
  {
    nodes.clear();
    cells.clear();
    faces.clear();
    infNodes.clear();
    boundaryNodes.clear();
    infFaces.clear();
    boundaryFaces.clear();
    faceCells.clear();
    convexHull.clear();
    neighborDomains.clear();
    sharedNodes.clear();
    sharedFaces.clear();
    points.clear();
  }

  //! Returns true if the tessellation is empty (not defined),
  //! false otherwise.
  virtual bool empty() const
  {
    return nodes.empty() and cells.empty() and faces.empty() and
       infNodes.empty() and boundaryNodes.empty() and infFaces.empty() and boundaryFaces.empty() and faceCells.empty() and
      convexHull.empty() and points.empty();
  }

  //! An array of (Dimension*numPoints) values containing components of
  //! generator positions. The components are stored in node-major order and
  //! the 0th component of the ith node appears in nodes[Dimension*i].
  std::vector<RealType> points;

  //! An array of (Dimension*numNodes) values containing components of
  //! node positions. The components are stored in node-major order and
  //! the 0th component of the ith node appears in nodes[Dimension*i].
  std::vector<RealType> nodes;

  //! This two-dimensional array defines the cell-face topology of the
  //! mesh. A cell has an arbitrary number of faces in 2D and 3D.
  //! cells[i][j] gives the index of the jth face of the ith cell.
  //! A negative face index indicates the actual face index is the 1's
  //! complement of the value (~cells[i][j]) and the face is oriented
  //! with an inward pointing normal for cells[i].
  std::vector<std::vector<int> > cells;

  //! This two-dimensional array defines the topology of the faces of the
  //! mesh. A face has an arbitrary number of nodes in 3D and 2 nodes in 2D.
  //! faces[i][j] gives the index of the jth node of the ith face.
  //! Nodes for a given face are arranged counterclockwise around the face
  //! viewed from the "positive" (outside) direction.
  std::vector<std::vector<unsigned> > faces;

  //! Indices of all nodes that are on the boundary of the tessellation.
  std::vector<unsigned> infNodes;

  //! Indices of all nodes that are on the boundary of the tessellation.
  std::vector<unsigned> boundaryNodes;

  //! Indices of all faces on the boundary of the tessellation.
  std::vector<unsigned> infFaces;

  //! Indices of all faces on the boundary of the tessellation.
  std::vector<unsigned> boundaryFaces;

  //! An array of cell indices for each face, i.e., the cells that share
  //! the face.
  //! For a given cell there will be either 1 or 2 cells -- the cases with 1
  //! cell indicate a face on a boundary of the tessellation.
  std::vector<std::vector<int> > faceCells;

  //! A PLC connecting the generating points belonging to the convex hull
  //! of the point distribution. Not all Tessellators hand back the convex
  //! hull, so this may be empty, in which case you must compute the convex
  //! hull yourself.
  PLC<Dimension> convexHull;

  //! Parallel data structure: the set of neighbor domains this portion of
  //! the tessellation is in contact with.
  std::vector<unsigned> neighborDomains;

  //! Parallel data structure: the nodes and faces this domain shares with
  //! each neighbor domain.
  //! NOTE: we implicitly assume that any domains of rank less than ours we
  //!       are receiving from, while any domains of greater rank we send
  //!       to.
  std::vector<std::vector<unsigned> > sharedNodes, sharedFaces;

  //! Find the set of cells that touch each mesh node.
  std::vector<std::set<unsigned> > computeNodeCells() {
    std::vector<std::set<unsigned> > result(nodes.size()/Dimension);
    for (auto i = 0u; i < cells.size(); ++i) {
      for (std::vector<int>::const_iterator faceItr = cells[i].begin();
           faceItr != cells[i].end();
           ++faceItr) {
        const unsigned iface = *faceItr < 0 ? ~(*faceItr) : *faceItr;
        for (std::vector<unsigned>::const_iterator nodeItr = faces[iface].begin();
             nodeItr != faces[iface].end();
             ++nodeItr) {
          POLY_ASSERT(*nodeItr < result.size());
          result[*nodeItr].insert(i);
        }
      }
    }
    return result;
  }


  //! Collect the nodes around each cell
  std::vector<std::set<unsigned> > computeCellToNodes()
  {
    std::vector<std::set<unsigned> > result(cells.size());
    for (unsigned i = 0; i != cells.size(); ++i){
      for (std::vector<int>::const_iterator faceItr = cells[i].begin();
           faceItr != cells[i].end(); ++faceItr){
        const unsigned iface = *faceItr < 0 ? ~(*faceItr) : *faceItr;
        POLY_ASSERT(iface < faceCells.size());
        for (std::vector<unsigned>::const_iterator nodeItr = faces[iface].begin();
             nodeItr != faces[iface].end(); ++nodeItr) {
          POLY_ASSERT(*nodeItr < nodes.size());
          result[i].insert(*nodeItr);
        }
      }
    }
    return result;
  }

  RealCell getCell(const unsigned cellIndex) const {
    return Cell<Dimension, RealType>::extractCell(nodes, cells[cellIndex], faces);
  }

  //! output operator.
  friend std::ostream& operator<<(std::ostream& s, const Tessellation& mesh) {
    for (int i = 0; i < mesh.cells.size(); ++i) {
      s << mesh.getCell(i);
    }
    return s;
  }

  // template<int D = Dimension>
  // std::enable_if_t<D == 2, void>
  // void computeCellCentroid(const unsigned ci,
  //                          RealType* ccent) const;

  template<int D = Dimension>
  std::enable_if_t<D == 2, void>
  computeCellCentroidAndSignedArea(const unsigned ci,
                                   const RealType& tol,
                                   RealType* ccent,
                                   RealType& area) const;
  template<int D = Dimension>
  std::enable_if_t<D == 3, void>
  computeFaceCentroidAndNormal(const unsigned fi,
                               RealType* fcent,
                               RealType* fhat) const;

  template<int D = Dimension>
  std::enable_if_t<D == 3, void>
  computeCellCentroidAndSignedVolume(const unsigned ci,
                                     RealType* ccent,
                                     RealType& cvol) const;
  void computeFaceCells() {
    faceCells.clear();
    auto numFaces = faces.size();
    auto numCells = cells.size();
    faceCells.resize(numFaces);
    // Build faceCells connectivity: for each cell, mark which faces touch it
    // Cells store signed face indices where negative means inverted orientation.
    for (unsigned i = 0; i < numCells; ++i) {
      const unsigned nf = cells[i].size();
      for (unsigned j = 0; j < nf; ++j) {
        auto k = cells[i][j];
        if (k < 0) {
          // Negative index: inverted face orientation
          POLY_ASSERT2(~k < numFaces, k << " " << ~k << " " << numFaces);
          faceCells[~k].push_back(~i);
        } else {
          // Positive index: normal face orientation
          POLY_ASSERT2(k < numFaces, k << " " << numFaces);
          faceCells[k].push_back(i);
        }
      }
    }
  }

private:

  // Disallowed.
  Tessellation(const Tessellation&);
  Tessellation& operator=(const Tessellation&);
};

}

#endif
