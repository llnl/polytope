// Try tessellating a simple lattice of generators in a box.

#include <iostream>
#include <vector>
#include <map>
#include <stdint.h>
#include <stdlib.h>

#include "polytope.hh"

#define POLY_CHECK(x) if (!(x)) { cout << "FAIL: " << #x << endl; exit(-1); }

using namespace std;

//------------------------------------------------------------------------------
// Hash a 3 dimensional position to a single 64 unsigned int.
//------------------------------------------------------------------------------
template<typename Real>
inline
uint64_t 
hash3position(const Real* pos,
              const Real* xmin,
              const Real* xmax,
              const Real& dx) {
  uint64_t result = static_cast<uint64_t>((max(xmin[2], min(xmax[2], pos[2])) - xmin[2])/dx + 0.5);
  result <<= 21;
  result |= static_cast<uint64_t>((max(xmin[1], min(xmax[1], pos[1])) - xmin[1])/dx + 0.5);
  result <<= 21;
  result |= static_cast<uint64_t>((max(xmin[0], min(xmax[0], pos[0])) - xmin[0])/dx + 0.5);
  return result;
}

//------------------------------------------------------------------------------
// The main test.
//------------------------------------------------------------------------------
int main() {
  const double x1 = 0.0, y1 = 0.0, z1 = 0.0;
  const double x2 = 100.0, y2 = 100.0, z2 = 100.0;

  // Try tessellating increasing numbers of generators.
  for (unsigned nx = 2; nx != 50; ++nx) {
    cout << "Testing nx=" << nx << endl;

    // Create the generators.
    vector<double> generators;
    const double dx = (x2 - x1)/nx, dy = (y2 - y1)/nx, dz = (z2 - z1)/nx;
    unsigned ix, iy, iz;
    double xi, yi, zi;
    for (iz = 0; iz != nx; ++iz) {
      zi = std::max(z1, std::min(z2, z1 + (iz + 0.5)*dz));
      for (iy = 0; iy != nx; ++iy) {
        yi = std::max(y1, std::min(y2, y1 + (iy + 0.5)*dy));
        for (ix = 0; ix != nx; ++ix) {
          xi = std::max(x1, std::min(x2, x1 + (ix + 0.5)*dx));
          generators.push_back(xi);
          generators.push_back(yi);
          generators.push_back(zi);
        }
      }
    }

    // Create the tessellation.
    double xmin[3] = { x1, y1, z1 };
    double xmax[3] = { x2, y2, z2 };
    polytope::Tessellation<3, double> mesh;
    polytope::VoroPP_3d<double> voro;
    voro.tessellate(generators, xmin, xmax, mesh);

    // Spew the mesh statistics.
    cout << "   num mesh nodes : " << mesh.nodes.size()/3 << endl;
    cout << "   num mesh cells : " << mesh.cells.size() << endl;
    cout << "   num mesh faces : " << mesh.faces.size() << endl;
//     cout << "Node positions: " << endl;
//     for (unsigned i = 0; i != mesh.nodes.size()/3; ++i) {
//       cout << "   Node " << i << " @ (" << mesh.nodes[3*i] << " " << mesh.nodes[3*i + 1] << " " << mesh.nodes[3*i + 2] << ")" << " " << hash3position(&mesh.nodes[3*i], xmin, xmax, 0.1*dx) 
//            << " " << static_cast<uint64_t>((max(xmin[1], min(xmax[1], mesh.nodes[3*i + 1])) - xmin[1])/dx + 0.5)
//            << " " << static_cast<uint64_t>((max(xmin[0], min(xmax[0], mesh.nodes[3*i + 0])) - xmin[0])/dx + 0.5)
//            << endl;

//     }
//     cout << "Face node sets: " << endl;
//     for (unsigned i = 0; i != nx*nx*nx; ++i) {
//       cout << "   FACES for mesh cell " << i << " :";
//       for (unsigned j = 0; j != mesh.cells[i].size(); ++j) cout << " " << mesh.cells[i][j];
//       cout << endl;
//     }
//     for (unsigned i = 0; i != mesh.faces.size(); ++i) {
//       double xf = 0.0, yf = 0.0, zf = 0.0;
//       cout << "   NODES for mesh face " << i << " :";
//       for (unsigned j = 0; j != mesh.faces[i].size(); ++j) {
//         unsigned k = mesh.faces[i][j];
//         cout << " " << k;
//         xf += mesh.nodes[3*k];
//         yf += mesh.nodes[3*k + 1];
//         zf += mesh.nodes[3*k + 2];
//       }
//       xf /= mesh.faces[i].size();
//       yf /= mesh.faces[i].size();
//       zf /= mesh.faces[i].size();
//       cout << " @ (" << xf << " " << yf << " " << zf << ")"  << endl;
//     }

    // Now do the checks.
    // If we bin the nodes up a small fraction of a cell size they should be unique!
    map<uint64_t, unsigned> nodeHash2ID;
    unsigned i;
    uint64_t hashi;
    for (i = 0; i != nx*nx*nx; ++i) {
      hashi = hash3position(&mesh.nodes[3*i], xmin, xmax, 0.1*dx);
      POLY_CHECK(nodeHash2ID.find(hashi) == nodeHash2ID.end());
      nodeHash2ID[hashi] = i;
    }

    // Check sizes.
    POLY_CHECK(mesh.nodes.size()/3 == (nx + 1)*(nx + 1)*(nx + 1));
    POLY_CHECK(mesh.cells.size() == nx*nx*nx);
    for (unsigned i = 0; i != nx*nx*nx; ++i) POLY_CHECK(mesh.cells[i].size() == 6);
    POLY_CHECK(mesh.faces.size() == 3*nx*nx*(nx + 1));
  }

  cout << "PASS" << endl;
  return 0;
}
