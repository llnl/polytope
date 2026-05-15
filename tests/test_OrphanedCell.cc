// test_OrphanedCell
//
// 3x3 Grid with Cartesian generators. Two rectangles jut out of the top
// and bottom boundaries, dividing the top-middle and bottom-middle cells.
// For each of the two cells: the divided part that contains the generator
// becomes the cell, while the other part becomes an orphaned cell.
//
// This tests the "cell adoption" capability for 2D tessellators
// which appropriates the area in the orphaned pieces to its neighboring
// cells while maintaining the Voronoi property locally
//
// The domain looks like this:
//   _____  ________
//  |     ||        |
//  |     ||        |
//  |     --        |
//  |               |
//  |     --        |
//  |     ||        |
//  |_____||________|


#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "polytope_test_utilities.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// test
// -----------------------------------------------------------------------
void test(Tessellator<2,double>& tessellator,
          const bool checkMesh) {

  // Initialize the bounded tessellation shtuff
  std::vector<double> PLCpoints;
  std::vector<double> points;
  PLC<2> boundary;
  Tessellation<2,double> mesh;

  // Test name for output
  string testName = "OrphanedCell_" + tessellator.name();
  
  // Boundary points
  PLCpoints.push_back(0.0);  PLCpoints.push_back(0.0);
  PLCpoints.push_back(1.2);  PLCpoints.push_back(0.0);
  PLCpoints.push_back(1.2);  PLCpoints.push_back(1.3);
  PLCpoints.push_back(1.3);  PLCpoints.push_back(1.3);
  PLCpoints.push_back(1.3);  PLCpoints.push_back(0.0);
  PLCpoints.push_back(3.0);  PLCpoints.push_back(0.0);
  PLCpoints.push_back(3.0);  PLCpoints.push_back(3.0);
  PLCpoints.push_back(1.3);  PLCpoints.push_back(3.0);
  PLCpoints.push_back(1.3);  PLCpoints.push_back(1.7);
  PLCpoints.push_back(1.2);  PLCpoints.push_back(1.7);
  PLCpoints.push_back(1.2);  PLCpoints.push_back(3.0);
  PLCpoints.push_back(0.0);  PLCpoints.push_back(3.0);
   
  // 3x3 Cartesian generators
  int ix, iy, nx = 3;
  double xi, yi;
  for (iy = 0; iy != nx; ++iy) {
    yi = iy + 0.5;
    for (ix = 0; ix != nx; ++ix) {
      xi = ix + 0.5;
      points.push_back(xi);  points.push_back(yi);
    }
  }
  
  // Hook together the facets
  int nSides = PLCpoints.size()/2;
  boundary.facets.resize( nSides, std::vector<int>(2) );
  for (unsigned i = 0; i != nSides; ++i){
    boundary.facets[i][0] = i;
    boundary.facets[i][1] = (i+1) % nSides;
  }
   
  // Tessellate
  tessellator.tessellate(points, PLCpoints, boundary, mesh);
  
  // Store the mesh using SiloWriter
  outputMesh(mesh, testName, points);

  //cerr << mesh << endl;

  if (checkMesh) {
    // Total mesh checks
    POLY_CHECK(mesh.cells.size()   == 9 );
    POLY_CHECK(mesh.nodes.size()/2 == 26);
    POLY_CHECK(mesh.faces.size()   == 34);
    // Individual cell checks
    POLY_CHECK(mesh.cells[0].size() == 5 );
    POLY_CHECK(mesh.cells[1].size() == 4 );
    POLY_CHECK(mesh.cells[2].size() == 4 );
    POLY_CHECK(mesh.cells[3].size() == 4 );
    POLY_CHECK(mesh.cells[4].size() == 12);
    POLY_CHECK(mesh.cells[5].size() == 4 );
    POLY_CHECK(mesh.cells[6].size() == 5 );
    POLY_CHECK(mesh.cells[7].size() == 4 );
    POLY_CHECK(mesh.cells[8].size() == 4 );
    // Tessellation area check
    const double trueArea = 8.74;
    const double tessArea = computeTessellationArea(mesh);
    const double fracerr  = std::abs(trueArea - tessArea)/trueArea;
    const double tol      = 1.0e-7;
    POLY_CHECK2(fracerr < tol, "Relative error in the tessellation "
                << "area exceeds tolerance:" << endl
                << "      Area = " << tessArea << endl
                << "     Error = " << trueArea - tessArea << endl
                << "Frac Error = " << fracerr);
  } else {
    cout << "Checking currently disabled for " << tessellator.name() << endl;
  }
}

// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#ifdef POLYTOPE_ENABLE_MPI
  MPI_Init(&argc, &argv);
#endif
   
#ifdef POLYTOPE_ENABLE_TRIANGLE
  {
    cout << "\nTriangle Tessellator:\n" << endl;
    TriangleTessellator<double> tessellator;
    test(tessellator, true);
  }
#endif

#ifdef POLYTOPE_ENABLE_BOOST
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator<double> tessellator;
    test(tessellator, true);
  }
#endif

  

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
