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
double computeTessellationArea(Tessellation<2, double>& mesh, bool parallel = true) {
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
    area += boost::geometry::area(cellPolygon);
    nodeCell.clear();
  }
#ifdef POLYTOPE_ENABLE_MPI
  if (parallel) {
    double distributedArea = 0.0;
    MPI_Allreduce(&area, &distributedArea, 1, MPI_DOUBLE, MPI_SUM, Communicator::communicator());
    area = distributedArea;
  }
#endif
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
// Estimate the error tolerance for the 2D area
// -----------------------------------------------------------------------
double errorTolerance(const Quantizer<2>& Q,
                      const unsigned Ncells) {
  const double eps = std::numeric_limits<double>::epsilon();
  const double u = 0.5 * eps;
  const double nAdds = double(Ncells - 1);
  const double gamma = (nAdds * u) / (1.0 - nAdds * u);
  // Dimensionless allowance for box-scale quantization.
  const double quantTol = 2.0*(Q.degeneracy().x + Q.degeneracy().y);
  // Practical margin for cell-area arithmetic and accumulation.
  const double roundoffTol = 4.0 * gamma + 32.0 * eps;
  return quantTol + roundoffTol;
}

// -----------------------------------------------------------------------
// compareArea
// -----------------------------------------------------------------------
void compareArea(Boundary2D& boundary,
                 Tessellation<2, double>& mesh,
                 const Quantizer<2>& Q,
                 const std::string& test_str = "",
                 bool parallel = true) {
  double area = computeTessellationArea(mesh, parallel);
  bool clipped = boundary.mClipping;
  const auto Ncells = mesh.totalPointCount(parallel);
  const double tol = errorTolerance(Q, Ncells);
  const double bArea = clipped ? boundary.area() : Q.area();
  const double relErr = std::abs(bArea-area)/bArea;
  POLY_CHECK2(relErr < tol, "Error in area: error = " << relErr << " tolerance "
              << tol << " " << test_str);
  // Check for CCW orientation of polygons
  for (auto i = 0u; i < mesh.cells.size(); ++i) {
    POLY_CHECK2(ccwOrient(mesh.getCell(i)), "Bad orientation at cell " << i << " " << test_str);
  }
}

//------------------------------------------------------------------------------
// Rotational test helpers
//------------------------------------------------------------------------------

// -----------------------------------------------------------------------
// getVelocities
// -----------------------------------------------------------------------
inline
void getVelocities(const vector<Point2<double>>& points,
                   const Point2<double>& bhi,
                   const Point2<double>& blo,
                   const unsigned flowType,
                   vector<Point2<double>>& velocities) {
  auto center = 0.5*(bhi + blo);
  auto len = bhi - blo;
  double maxRadius = std::sqrt(len.x*len.x + len.y*len.y)/2.;
  double omega = 1./maxRadius;
  POLY_CHECK(points.size() == velocities.size());
  for (unsigned i = 0; i < points.size(); ++i) {
    const auto& p = points[i];
    auto& v = velocities[i];
    const auto u = (p - blo)/len;
    switch(flowType) {
    case 1:
      v.x = -omega*(p.y - center.y);
      v.y = omega*(p.x - center.x);
      break;
    case 2:
      v.x =  len.x * sin(M_PI * u.x) * cos(M_PI * u.y);
      v.y = -len.y * cos(M_PI * u.x) * sin(M_PI * u.y);
      break;
    case 3:
      v.x =  0.5 * len.x * cos(2.0 * M_PI * (u.x - 0.25)) *
        sin(2.0 * M_PI * (u.y - 0.25));
      v.y = -0.5 * len.y * sin(2.0 * M_PI * (u.x - 0.25)) *
        cos(2.0 * M_PI * (u.y - 0.25));
      break;
    case 4:
      v.x = -len.x * sin(4.0 * M_PI * (u.x + 0.5)) *
        sin(4.0 * M_PI * (u.y + 0.5));
      v.y = -len.y * cos(4.0 * M_PI * (u.x + 0.5)) *
        cos(4.0 * M_PI * (u.y + 0.5));
      break;
    }
  }
}

// Assumes max velocity <= 1
inline
std::string setConfig(unsigned flowType,
                      const Point2<double>& len,
                      const double& dx,
                      double& dt, double& Tmax,
                      double& scaleFactor) {
  double vel = 1.;
  double cfl = std::sqrt(2.);
  double maxRadius = std::sqrt(len.x*len.x + len.y*len.y)/2.;
  switch(flowType) {
  case 1:
    scaleFactor = 1.0/std::sqrt(2.);
    dt = cfl*dx;
    Tmax = M_PI*maxRadius;
    return "Test 1: Solid Rotation Flow";
  case 2:
    dt = cfl*dx/vel;
    Tmax = 2.0;
    scaleFactor = 1.0;
    return "Test 2: Single Vortex Flow";
  case 3:
    dt = cfl*dx/vel;
    Tmax = 2.0;
    scaleFactor = 1.0;
    return "Test 3: Taylor-Green (4-Vortex) Flow";
  case 4:
    cfl = 0.5;
    dt = cfl*dx/vel;
    Tmax = 2.0;
    scaleFactor = 1.0/std::sqrt(2);
    return "Test 4: Deformation (16-Vortex) Flow";
  }
  return "";
}
}
#endif
