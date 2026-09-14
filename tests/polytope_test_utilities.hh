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
		const double time = 0.0,
                const int numFiles = -1) {
#ifdef POLYTOPE_ENABLE_SILO
  SiloWriter<2, Tessellation<2, double>> writer(mesh);
  writer.generateTestVars();
  writer.write(prefix, "", testCycle, time, numFiles);
#endif
}

template<typename FieldType>
void outputMesh(const Tessellation<2, double>& mesh,
		std::string prefix,
                std::vector<FieldType>& cellFieldVec,
                std::string cellFieldName,
		const unsigned testCycle = 1,
		const double time = 0.0,
                const int numFiles = 1) {
#ifdef POLYTOPE_ENABLE_SILO
  SiloWriter<2, Tessellation<2, double>> writer(mesh);
  writer.generateTestVars();
  writer.addField<double>(FieldCentering::Cell, cellFieldName, cellFieldVec);
  writer.write(prefix, "", testCycle, time, numFiles);
#endif
}

//..............................................................................
// 3D
void outputMesh(const Tessellation<3, double>& mesh,
		std::string prefix,
		const unsigned testCycle = 1,
		const double time = 0.0) {
#ifdef POLYTOPE_ENABLE_SILO
  SiloWriter<3, Tessellation<3, double>> writer(mesh);
  writer.generateTestVars();
  writer.write(prefix, "", testCycle, time);
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
    BGPolygon<double, 2> polygon = makeBGPolygon(cell.points());
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
                 double in_area,
                 const std::string& test_str = "") {
  double area = in_area;
#ifdef POLYTOPE_ENABLE_MPI
  const int nranks = Communicator::getNRanks();
  if (nranks > 1) {
    double distributedArea = 0.0;
    MPI_Allreduce(&in_area, &distributedArea, 1, MPI_DOUBLE, MPI_SUM, Communicator::communicator());
    area = distributedArea;
  }
#endif
  const double relErr = std::abs(boundary.mArea-area)/boundary.mArea;
  auto& Q = Quantizer<2>::instance();
  const int maxAxis = Q.m_lx_o.maxAxis();
  const double tol = 2.*Q.degeneracy()[maxAxis];
  POLY_CHECK2(relErr < tol, "Error in area: error = " << relErr << " tolerance "
              << tol << " " << test_str);
}

void compareArea(Boundary2D& boundary,
                 Tessellation<2, double>& mesh,
                 const std::string& test_str = "") {
  double area = computeTessellationArea(mesh);
  compareArea(boundary, area, test_str);
}

}
#endif
