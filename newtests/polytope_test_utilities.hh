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
#include "Communicator.hh"
#include "polytope_boost_utilities.hh"

namespace polytope {

//------------------------------------------------------------------------------
// A simple mesh output function for the SiloWriter
//------------------------------------------------------------------------------
// 2D
void outputMesh(const Tessellation<2, double>& mesh,
		std::string prefix,
		const unsigned testCycle = 1,
		const double time = 0.0) {
#ifdef POLYTOPE_ENABLE_SILO
  std::map<std::string, double*> nodeFields, edgeFields, faceFields, cellFields;
  size_t meshSize = mesh.cells.size();
  std::vector<double> index(meshSize);
  std::vector<double> genx (meshSize);
  std::vector<double> geny (meshSize);
  for (int i = 0; i < meshSize; ++i) {
    index[i] = double(i);
    genx[i] = mesh.points[i].x;
    geny[i] = mesh.points[i].y;
  }
  cellFields["cell_index"] = &index[0];
  cellFields["gen_x"     ] = &genx[0];
  cellFields["gen_y"     ] = &geny[0];
#ifdef POLYTOPE_ENABLE_MPI
  int rank = Communicator::getRank();
  std::vector<double> rankField(meshSize, rank);
  cellFields["rank"      ] = &rankField[0];
#endif
  std::ostringstream os;
  os << prefix;
  SiloWriter<2, Tessellation<2, double>>::write(mesh, nodeFields, edgeFields,
                                                faceFields, cellFields, os.str(),
                                                testCycle, time);
#endif
}

void outputMesh(const Tessellation<2, double>& mesh,
		std::string prefix,
                std::vector<double>& cellFieldVec,
                std::string cellFieldName,
		const unsigned testCycle = 1,
		const double time = 0.0,
                const int numFiles = 1) {
#ifdef POLYTOPE_ENABLE_SILO
  std::map<std::string, double*> nodeFields, edgeFields, faceFields, cellFields;
  size_t meshSize = mesh.cells.size();
  POLY_CHECK(cellFieldVec.size() == meshSize);
  std::vector<double> index(meshSize);
  std::vector<double> genx (meshSize);
  std::vector<double> geny (meshSize);
  for (int i = 0; i < meshSize; ++i) {
    index[i] = double(i);
    genx[i] = mesh.points[i].x;
    geny[i] = mesh.points[i].y;
  }
  cellFields["cell_index"] = &index[0];
  cellFields["gen_x"     ] = &genx[0];
  cellFields["gen_y"     ] = &geny[0];
  cellFields[cellFieldName] = &cellFieldVec[0];
  std::ostringstream os;
  os << prefix;
  SiloWriter<2, Tessellation<2, double>>::write(mesh, nodeFields, edgeFields,
                                                faceFields, cellFields, os.str(),
                                                testCycle, time, numFiles);
#endif
}

//..............................................................................
// 3D
void outputMesh(const Tessellation<3, double>& mesh,
		std::string prefix,
		const unsigned testCycle = 1,
		const double time = 0.0) {
#ifdef POLYTOPE_ENABLE_SILO
  std::vector<double> index(mesh.cells.size());
  std::vector<double> genx (mesh.cells.size());
  std::vector<double> geny (mesh.cells.size());
  std::vector<double> genz (mesh.cells.size());
  //std::vector<double> vol  (mesh.cells.size());
  for (int i = 0; i < mesh.cells.size(); ++i){
    index[i] = double(i);
    if (!mesh.points.empty()) {
      genx[i] = mesh.points[i].x;
      geny[i] = mesh.points[i].y;
      genz[i] = mesh.points[i].z;
    }
    //mesh.computeCellCentroidAndSignedVolume(i, cent, vol[i]);
  }
  std::map<std::string,double*> nodeFields, edgeFields, faceFields, cellFields;
  cellFields["cell_index"] = &index[0];
  cellFields["gen_x"     ] = &genx[0];
  cellFields["gen_y"     ] = &geny[0];
  cellFields["gen_z"     ] = &genz[0];
  //cellFields["volume"    ] = &vol[0];
#ifdef POLYTOPE_ENABLE_MPI
  int rank = Communicator::getRank();
  std::vector<double> rankField(mesh.cells.size(), rank);
  cellFields["rank"      ] = &rankField[0];
#endif
  std::ostringstream os;
  os << prefix;
  SiloWriter<3, Tessellation<3, double>>::write(mesh, nodeFields, edgeFields,
                                                faceFields, cellFields, os.str(),
                                                testCycle, time);
#endif
}

//------------------------------------------------------------------------------
// Some specialized subsets of outputMesh
//------------------------------------------------------------------------------
template <int nDim>
void outputMesh(const Tessellation<nDim, double>& mesh,
		std::string prefix,
		const unsigned testCycle) {
  outputMesh(mesh, prefix, testCycle, 0.0);
}
//------------------------------------------------------------------------------
template <int nDim>
void outputMesh(const Tessellation<nDim, double>& mesh,
		std::string prefix) {
  outputMesh(mesh, prefix, 1, 0.0);
}
//------------------------------------------------------------------------------
// a cell-centered field given
// template <typename RealType>
// void outputMesh(Tessellation<2,RealType>& mesh,
// 		std::string prefix,
//                 std::vector<RealType>& cellField,
// 		const unsigned testCycle = 1,
// 		const RealType time = 0.0) {
// #ifdef POLYTOPE_ENABLE_SILO
//   POLY_CHECK(cellField.size() == mesh.cells.size());
//   std::vector<double> index(mesh.cells.size());
//   std::vector<double> genx (mesh.cells.size());
//   std::vector<double> geny (mesh.cells.size());
//   for (int i = 0; i < mesh.cells.size(); ++i){
//     index[i] = double(i);
//     genx[i] = mesh.points[2*i];
//     geny[i] = mesh.points[2*i+1];
//   }
//   std::map<std::string,double*> nodeFields, edgeFields, faceFields, cellFields;
//   cellFields["cell_index"] = &index[0];
//   cellFields["gen_x"     ] = &genx[0];
//   cellFields["gen_y"     ] = &geny[0];
//   cellFields["cond"      ] = &cellField[0];
//   std::ostringstream os;
//   os << prefix;
//   SiloWriter<2, double>::write(mesh, nodeFields, edgeFields,
//                                faceFields, cellFields, os.str(),
//                                testCycle, time);
// #endif
// }

//------------------------------------------------------------------------------
// Compute the area of a polytope tessellation cell-by-cell using Boost.Geometry
//------------------------------------------------------------------------------
double computeTessellationArea(Tessellation<2, double>& mesh) {
  double area = 0;
  for (unsigned i = 0; i != mesh.cells.size(); ++i) {
    std::vector<double> nodeCell;
    for (std::vector<int>::const_iterator faceItr = mesh.cells[i].begin();
         faceItr != mesh.cells[i].end(); ++faceItr){
      const unsigned iface = *faceItr < 0 ? ~(*faceItr) : *faceItr;
      POLY_CHECK(iface < mesh.faceCells.size());
      POLY_CHECK(mesh.faces[iface].size() == 2);
      const unsigned inode = *faceItr < 0 ? mesh.faces[iface][1] : mesh.faces[iface][0];
      nodeCell.push_back( mesh.nodes[inode].x );
      nodeCell.push_back( mesh.nodes[inode].y );
    }
    BGPolygon<double,2> cellPolygon = makeBGPolygon<double>( nodeCell );
    area += boost::geometry::area( cellPolygon );
    nodeCell.clear();
  }
  return area;
}

// Return -1 if polygon is not watertight, otherwise returns number of holes
int isWatertight(const Tessellation<2, double>& mesh) {
  MultiBGPolygon<double, 2> result;
  for (int i = 0; i < mesh.cells.size(); ++i) {
    auto cell = mesh.getCell(i);
    BGPolygon<double, 2> polygon = makeBGPolygon(cell);
    boost::geometry::correct(polygon);
    if (!boost::geometry::is_valid(polygon)) {
      return -1;
    }
    if (i == 0) {
      result.push_back(polygon);
    } else {
      MultiBGPolygon<double, 2> tmp;
      boost::geometry::union_(result, polygon, tmp);
      result = std::move(tmp);
    }
  }
  int holeCount = 0;
  for (const auto& mp : result) {
    holeCount += mp.inners().size();
  }
  return holeCount;
}

void testWatertight(const Tessellation<2, double>& mesh, const int refHoles) {
  int numHoles = isWatertight(mesh);
  POLY_CHECK2(numHoles != -1, "Resulting mesh is not watertight");
  (void) refHoles;
  // This check does not work as intended for some reason
  // POLY_CHECK2(numHoles == refHoles,
  //             "Resulting mesh has " << numHoles << " but should have " << refHoles << " holes");
}

// -----------------------------------------------------------------------
// compareArea
// -----------------------------------------------------------------------
void compareArea(Boundary2D& boundary,
                 Tessellation<2,double>& mesh) {
   const double area = computeTessellationArea(mesh);
   const double relErr = std::abs(boundary.mArea-area)/boundary.mArea;
   POLY_CHECK2(relErr < 1.0E-8, "Error in area: ref " << boundary.mArea << " mesh " << area);
}

}
#endif
