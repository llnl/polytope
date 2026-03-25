// Try tessellating a simple lattice of generators in a box.

#include <iostream>
#include <vector>
#include <cassert>
#include <ctime>

#include "polytope.hh"
#include "polytope_test_utilities.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

//------------------------------------------------------------------------------
// Emergency dump of the mesh.
//------------------------------------------------------------------------------
std::string
escapePod(const unsigned nx,
          const vector<double>& generators,
          const Tessellation<3, double>& mesh) {
    std::stringstream os;
    os << "test_TetgenTessellator_" << nx << "x" << nx << "x" << nx;
    outputMesh(mesh, os.str(), generators);
    return " : attempted to write to file " + os.str();
}

//------------------------------------------------------------------------------
// unbounded.
//------------------------------------------------------------------------------
void unboundedTessellation(const unsigned nx,
                           const vector<double>& generators) {

  // Create the tessellation.
  Tessellation<3, double> mesh;
  TetgenTessellator tetgen;
  tetgen.tessellate(generators, mesh);

  // escapePod(nx, generators, mesh);

  // Check for validity.
  const unsigned nx1 = nx - 1;
  POLY_CHECK(mesh.nodes.size()/3 == nx1*nx1*nx1 + 6*nx1*nx1);
  POLY_CHECK(mesh.cells.size() == nx*nx*nx);
  POLY_CHECK2(mesh.infNodes.size() == 6*nx1*nx1, "Number of infNodes: " << mesh.infNodes.size() << escapePod(nx, generators, mesh));
  for (unsigned i = 0; i != nx*nx*nx; ++i) {
    const unsigned 
      ix = i % nx,
      iy = (i / nx) % nx,
      iz = i / (nx*nx);
    const unsigned ntouch = (unsigned(ix == 0 or ix == nx1) + 
                             unsigned(iy == 0 or iy == nx1) + 
                             unsigned(iz == 0 or iz == nx1));
    if (ntouch == 3) {
      // Corner cell.
      POLY_CHECK2(mesh.cells[i].size() == 4, escapePod(nx, generators, mesh));
    } else if (ntouch == 2) {
      // Along one of the edges of the volume.
      POLY_CHECK2(mesh.cells[i].size() == 5, escapePod(nx, generators, mesh));
    } else if (ntouch == 1) {
      // Along one of the faces of the volume.
      POLY_CHECK2(mesh.cells[i].size() == 6, escapePod(nx, generators, mesh));
    } else {
      // Interior, fully bounded cell.
      POLY_CHECK2(mesh.cells[i].size() == 6, escapePod(nx, generators, mesh));
    }

    // Check face orientations.
    for (unsigned j = 0; j != mesh.cells[i].size(); ++j) {
      const int faceSgn = geometry::sgn(mesh.cells[i][j]);
      const unsigned iface = internal::positiveID(mesh.cells[i][j]);
      const unsigned nnodes = mesh.faces[iface].size();
      for (unsigned k = 0; k != nnodes; ++k) {
        const unsigned inode1 = mesh.faces[iface][k];
        const unsigned inode2 = mesh.faces[iface][(k + 1) % nnodes];
        const unsigned inode3 = mesh.faces[iface][(k + 2) % nnodes];
        POLY_CHECK(inode1 < mesh.nodes.size()/3);
        POLY_CHECK(inode2 < mesh.nodes.size()/3);
        POLY_CHECK(inode3 < mesh.nodes.size()/3);
        const double vol = geometry::tetrahedralVolume6(&generators[3*i],
                                                        &mesh.nodes[3*inode3],
                                                        &mesh.nodes[3*inode2],
                                                        &mesh.nodes[3*inode1]);
        POLY_CHECK2(faceSgn*vol > 0.0,
                    "Volume sign error : " << i << " " << iface << " "
                    << faceSgn << " " << vol << " : " << nnodes << " " << inode1 << " " << inode2 << " " << inode3 << " : "
                    << " (" << generators[3*i] << " " << generators[3*i+1] << " " << generators[3*i+2] << ") "
                    << " (" << mesh.nodes[3*inode1] << " " << mesh.nodes[3*inode1+1] << " " << mesh.nodes[3*inode1+2] << ") "
                    << " (" << mesh.nodes[3*inode2] << " " << mesh.nodes[3*inode2+1] << " " << mesh.nodes[3*inode2+2] << ") "
                    << " (" << mesh.nodes[3*inode3] << " " << mesh.nodes[3*inode3+1] << " " << mesh.nodes[3*inode3+2] << ") ");
      }
    }
  }

}

//------------------------------------------------------------------------------
// bounded by a box.
//------------------------------------------------------------------------------
void boxBoundedTessellation(const unsigned nx,
                            const double x1, const double y1, const double z1,
                            const double x2, const double y2, const double z2,
                            const vector<double>& generators) {

  // Create the tessellation.
  Tessellation<3, double> mesh;
  TetgenTessellator tetgen;
  double low[3]  = {x1, y1, z1};
  double high[3] = {x2, y2, z2};
  tetgen.tessellate(generators, low, high, mesh);
  escapePod(nx, generators, mesh);

  // Check for validity.
  const unsigned nx1 = nx + 1;
  // POLY_CHECK(mesh.nodes.size()/3 == nx1*nx1*nx1);
  POLY_CHECK(mesh.cells.size() == nx*nx*nx);
  POLY_CHECK(mesh.infNodes.size() == 0);
  POLY_CHECK(mesh.infFaces.size() == 0);

  // Check nodes.
  const double tol = 1.0e-5*(x2 - x1);
  for (unsigned i = 0; i != mesh.nodes.size()/3; ++i) {
    POLY_CHECK2(x1 - tol <= mesh.nodes[3*i  ] and mesh.nodes[3*i  ] <= x2 + tol, "x bounds error: " << mesh.nodes[3*i  ] << " !\\in [" << x1 << ", " << x2 << "]");
    POLY_CHECK2(y1 - tol <= mesh.nodes[3*i+1] and mesh.nodes[3*i+1] <= y2 + tol, "y bounds error: " << mesh.nodes[3*i+1] << " !\\in [" << y1 << ", " << y2 << "]");
    POLY_CHECK2(z1 - tol <= mesh.nodes[3*i+2] and mesh.nodes[3*i+2] <= z2 + tol, "z bounds error: " << mesh.nodes[3*i+2] << " !\\in [" << z1 << ", " << z2 << "]");
  }

  // Check cells.
  for (unsigned i = 0; i != nx*nx*nx; ++i) {
    POLY_CHECK2(mesh.cells[i].size() == 6, escapePod(nx, generators, mesh));

    // Check face orientations.
    for (unsigned j = 0; j != mesh.cells[i].size(); ++j) {
      const int faceSgn = geometry::sgn(mesh.cells[i][j]);
      const unsigned iface = internal::positiveID(mesh.cells[i][j]);
      const unsigned nnodes = mesh.faces[iface].size();
      for (unsigned k = 0; k != nnodes; ++k) {
        const unsigned inode1 = mesh.faces[iface][k];
        const unsigned inode2 = mesh.faces[iface][(k + 1) % nnodes];
        const unsigned inode3 = mesh.faces[iface][(k + 2) % nnodes];
        POLY_CHECK(inode1 < mesh.nodes.size()/3);
        POLY_CHECK(inode2 < mesh.nodes.size()/3);
        POLY_CHECK(inode3 < mesh.nodes.size()/3);
        const double vol = geometry::tetrahedralVolume6(&generators[3*i],
                                                        &mesh.nodes[3*inode3],
                                                        &mesh.nodes[3*inode2],
                                                        &mesh.nodes[3*inode1]);
        POLY_CHECK2(faceSgn*vol > 0.0,
                    "Volume sign error : " << i << " " << iface << " "
                    << faceSgn << " " << vol << " : " << nnodes << " " << inode1 << " " << inode2 << " " << inode3 << " : "
                    << " (" << generators[3*i] << " " << generators[3*i+1] << " " << generators[3*i+2] << ") "
                    << " (" << mesh.nodes[3*inode1] << " " << mesh.nodes[3*inode1+1] << " " << mesh.nodes[3*inode1+2] << ") "
                    << " (" << mesh.nodes[3*inode2] << " " << mesh.nodes[3*inode2+1] << " " << mesh.nodes[3*inode2+2] << ") "
                    << " (" << mesh.nodes[3*inode3] << " " << mesh.nodes[3*inode3+1] << " " << mesh.nodes[3*inode3+2] << ") ");
      }
    }
  }

}

//------------------------------------------------------------------------------
// bounded by a PLC.
//------------------------------------------------------------------------------
void plcBoundedTessellation(const vector<double>& generators) {

  // Create the piecewise linear complex representing the box. Note that 
  // the box consists of facets that are defined by their connections to 
  // generating points.
  // Should look like the following:
  //
  //        6--------7            y
  //       /        /|            |
  //      /        / |            |
  //     2--------3  |             ------x
  //     |  .     |  |           /
  //     |  4.....|..5          z
  //     | .      | / 
  //     |.       |/
  //     0--------1             
  //
  // Create the vertices for our bounding surface.
  // vector<double> PLCpoints(3*8);
  // PLCpoints[3*0+0] = x1; PLCpoints[3*0+1] = y1; PLCpoints[3*0+2] = z2;
  // PLCpoints[3*1+0] = x2; PLCpoints[3*1+1] = y1; PLCpoints[3*1+2] = z2;
  // PLCpoints[3*2+0] = x1; PLCpoints[3*2+1] = y2; PLCpoints[3*2+2] = z2;
  // PLCpoints[3*3+0] = x2; PLCpoints[3*3+1] = y2; PLCpoints[3*3+2] = z2;
  // PLCpoints[3*4+0] = x1; PLCpoints[3*4+1] = y1; PLCpoints[3*4+2] = z1;
  // PLCpoints[3*5+0] = x2; PLCpoints[3*5+1] = y1; PLCpoints[3*5+2] = z1;
  // PLCpoints[3*6+0] = x1; PLCpoints[3*6+1] = y2; PLCpoints[3*6+2] = z1;
  // PLCpoints[3*7+0] = x2; PLCpoints[3*7+1] = y2; PLCpoints[3*7+2] = z1;

  // // 6 facets
  // PLC<3, double> box;
  // box.facets.resize(6);

  // // facet 0 -- bottom face.
  // box.facets[0].resize(4);
  // box.facets[0][0] = 0;
  // box.facets[0][1] = 4;
  // box.facets[0][2] = 5;
  // box.facets[0][3] = 1;

  // // facet 1 -- top face.
  // box.facets[1].resize(4);
  // box.facets[1][0] = 2;
  // box.facets[1][1] = 3;
  // box.facets[1][2] = 7;
  // box.facets[1][3] = 6;

  // // facet 2 -- left face.
  // box.facets[2].resize(4);
  // box.facets[2][0] = 0;
  // box.facets[2][1] = 2;
  // box.facets[2][2] = 6;
  // box.facets[2][3] = 4;

  // // facet 3 -- right face.
  // box.facets[3].resize(4);
  // box.facets[3][0] = 1;
  // box.facets[3][1] = 5;
  // box.facets[3][2] = 7;
  // box.facets[3][3] = 3;

  // // facet 4 -- front face.
  // box.facets[4].resize(4);
  // box.facets[4][0] = 0;
  // box.facets[4][1] = 1;
  // box.facets[4][2] = 3;
  // box.facets[4][3] = 2;

  // // facet 5 -- back face.
  // box.facets[5].resize(4);
  // box.facets[5][0] = 5;
  // box.facets[5][1] = 4;
  // box.facets[5][2] = 6;
  // box.facets[5][3] = 7;

}

// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv) {

#ifdef POLYTOPE_ENABLE_MPI
   MPI_Init(&argc, &argv);
#endif

  // Create the generators.
  const double x1 = 0.0, y1 = 0.0, z1 = 0.0;
  const double x2 = 1.0, y2 = 1.0, z2 = 1.0;
  // const double x1 = 0.0, y1 = 0.0, z1 = 0.0;
  // const double x2 = 100.0, y2 = 100.0, z2 = 100.0;
  unsigned ix, iy, iz;
  double xi, yi, zi;
  for (int nx = 2; nx != 30; ++nx) {
    cout << "============================== nx = " << nx << " ==============================" << endl;
    vector<double> generators;
    const double dx = (x2 - x1)/nx, dy = (y2 - y1)/nx, dz = (z2 - z1)/nx;
    for (iz = 0; iz != nx; ++iz) {
      zi = z1 + (iz + 0.5)*dz;
      for (iy = 0; iy != nx; ++iy) {
        yi = y1 + (iy + 0.5)*dy;
        for (ix = 0; ix != nx; ++ix) {
          xi = x1 + (ix + 0.5)*dx;
          generators.push_back(xi);
          generators.push_back(yi);
          generators.push_back(zi);
        }
      }
    }

    // Unbounded test.
    {
      cout << "Unbounded tessellation...";
      clock_t t0 = clock();
      unboundedTessellation(nx, generators);
      clock_t t1 = clock();
      cout << "required " << double(t1 - t0)/CLOCKS_PER_SEC << " seconds." << endl;
    }

    // Box bounded test.
    {
      cout << "Box bounded tessellation...";
      clock_t t0 = clock();
      boxBoundedTessellation(nx, x1, y1, z1, x2, y2, z2, generators);
      clock_t t1 = clock();
      cout << "required " << double(t1 - t0)/CLOCKS_PER_SEC << " seconds." << endl;
    }
  }

  

#ifdef POLYTOPE_ENABLE_MPI
   MPI_Finalize();
#endif
  return 0;
}
