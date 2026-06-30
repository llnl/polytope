//------------------------------------------------------------------------------
// A collection of random stuff useful for testing in polytope.
//------------------------------------------------------------------------------
#ifndef __polytope_test_utilities__
#define __polytope_test_utilities__

#include <sstream>
#include "polytope.hh"

#include "SiloWriter.hh"
#include "Tessellation.hh"
#include "Boundary2D.hh"
#include "Generators.hh"
#include "Tessellator.hh"

namespace polytope {

//------------------------------------------------------------------------------
// A simple mesh output function for the SiloWriter
//------------------------------------------------------------------------------
// 2D
template <typename RealType>
void outputMesh(const Tessellation<2,RealType>& mesh,
		std::string prefix,
		const unsigned testCycle = 1,
		const RealType time = 0.0) {
#ifdef POLYTOPE_ENABLE_SILO
  std::vector<double> index(mesh.cells.size());
  std::vector<double> genx (mesh.cells.size());
  std::vector<double> geny (mesh.cells.size());
  for (int i = 0; i < mesh.cells.size(); ++i) {
    index[i] = double(i);
    genx[i] = mesh.points[2*i];
    geny[i] = mesh.points[2*i+1];
  }
  std::map<std::string,double*> nodeFields, edgeFields, faceFields, cellFields;
  cellFields["cell_index"] = &index[0];
  cellFields["gen_x"     ] = &genx[0];
  cellFields["gen_y"     ] = &geny[0];
  std::ostringstream os;
  os << prefix;
  SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
                               faceFields, cellFields, os.str(),
                               testCycle, time);
#endif
}

//..............................................................................
// 3D
template <typename RealType>
void outputMesh(const Tessellation<3,RealType>& mesh,
		std::string prefix,
		const unsigned testCycle = 1,
		const RealType time = 0.0) {
#ifdef POLYTOPE_ENABLE_SILO
  std::vector<double> index(mesh.cells.size());
  std::vector<double> genx (mesh.cells.size());
  std::vector<double> geny (mesh.cells.size());
  std::vector<double> genz (mesh.cells.size());
  //std::vector<double> vol  (mesh.cells.size());
  double cent[3];
  for (int i = 0; i < mesh.cells.size(); ++i){
    index[i] = double(i);
    if (!mesh.points.empty()) {
      genx[i] = mesh.points[3*i  ];
      geny[i] = mesh.points[3*i+1];
      genz[i] = mesh.points[3*i+2];
    }
    //mesh.computeCellCentroidAndSignedVolume(i, cent, vol[i]);
  }
  std::map<std::string,double*> nodeFields, edgeFields, faceFields, cellFields;
  cellFields["cell_index"] = &index[0];
  cellFields["gen_x"     ] = &genx[0];
  cellFields["gen_y"     ] = &geny[0];
  cellFields["gen_z"     ] = &genz[0];
  //cellFields["volume"    ] = &vol[0];
  std::ostringstream os;
  os << prefix;
  SiloWriter<3, double>::write(mesh, nodeFields, edgeFields, 
                               faceFields, cellFields, os.str(),
                               testCycle, time);
#endif
}

//------------------------------------------------------------------------------
// Some specialized subsets of outputMesh
//------------------------------------------------------------------------------
template <int nDim, typename RealType>
void outputMesh(const Tessellation<nDim,RealType>& mesh,
		std::string prefix,
		const unsigned testCycle) {
  outputMesh(mesh, prefix, testCycle, 0.0);
}
//------------------------------------------------------------------------------
template <int nDim, typename RealType>
void outputMesh(const Tessellation<nDim,RealType>& mesh,
		std::string prefix) {
  outputMesh(mesh, prefix, 1, 0.0);
}
//------------------------------------------------------------------------------
// a cell-centered field given
template <typename RealType>
void outputMesh(Tessellation<2,RealType>& mesh,
		std::string prefix,
                std::vector<RealType>& cellField,
		const unsigned testCycle = 1,
		const RealType time = 0.0) {
#ifdef POLYTOPE_ENABLE_SILO
  POLY_CHECK(cellField.size() == mesh.cells.size());
  std::vector<double> index(mesh.cells.size());
  std::vector<double> genx (mesh.cells.size());
  std::vector<double> geny (mesh.cells.size());
  for (int i = 0; i < mesh.cells.size(); ++i){
    index[i] = double(i);
    genx[i] = mesh.points[2*i];
    geny[i] = mesh.points[2*i+1];
  }
  std::map<std::string,double*> nodeFields, edgeFields, faceFields, cellFields;
  cellFields["cell_index"] = &index[0];
  cellFields["gen_x"     ] = &genx[0];
  cellFields["gen_y"     ] = &geny[0];
  cellFields["cond"      ] = &cellField[0];
  std::ostringstream os;
  os << prefix;
  SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
                               faceFields, cellFields, os.str(),
                               testCycle, time);
#endif
}

//------------------------------------------------------------------------------
// Compute the area of a polytope tessellation cell-by-cell using Boost.Geometry
//------------------------------------------------------------------------------
template<typename RealType>
RealType computeTessellationArea(Tessellation<2,RealType>& mesh) {
  RealType area = 0;
  for (unsigned i = 0; i != mesh.cells.size(); ++i) {
    std::vector<RealType> nodeCell;
    for (std::vector<int>::const_iterator faceItr = mesh.cells[i].begin();
         faceItr != mesh.cells[i].end(); ++faceItr){
      const unsigned iface = *faceItr < 0 ? ~(*faceItr) : *faceItr;
      POLY_CHECK(iface < mesh.faceCells.size());
      POLY_CHECK(mesh.faces[iface].size() == 2);
      const unsigned inode = *faceItr < 0 ? mesh.faces[iface][1] : mesh.faces[iface][0];
      nodeCell.push_back( mesh.nodes[2*inode  ] );
      nodeCell.push_back( mesh.nodes[2*inode+1] );
    }
    BGPolygon<RealType,2> cellPolygon = makeBGPolygon<RealType>( nodeCell );
    area += boost::geometry::area( cellPolygon );
    nodeCell.clear();
  }
  return area;
}

}

#endif
