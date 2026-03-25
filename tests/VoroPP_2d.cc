// Try tessellating a simple lattice of generators in a box.

#include <iostream>
#include <vector>
#include "polytope.hh"
#include <stdlib.h>

#define POLY_CHECK(x) if (!(x)) { cout << "FAIL: " << #x << endl; exit(-1); }

using namespace std;

int main() {
  const double x1 = 0.0, y1 = 0.0;
  const double x2 = 1.0, y2 = 1.0;

  // Try tessellating increasing numbers of generators.
  for (unsigned nx = 2; nx != 100; ++nx) {
    cout << "Testing nx=" << nx << endl;

    // Create the generators.
    vector<double> generators;
    const double dx = (x2 - x1)/nx, dy = (y2 - y1)/nx;
    unsigned ix, iy;
    double xi, yi;
    for (iy = 0; iy != nx; ++iy) {
      yi = std::max(y1, std::min(y2, y1 + (iy + 0.5)*dy));
      for (ix = 0; ix != nx; ++ix) {
        xi = std::max(x1, std::min(x2, x1 + (ix + 0.5)*dx));
        generators.push_back(xi);
        generators.push_back(yi);
      }
    }

    // Create the tessellation.
    double xmin[2] = { x1, y1 };
    double xmax[2] = { x2, y2 };
    polytope::Tessellation<2, double> mesh;
    polytope::VoroPP_2d<double> voro;
    voro.tessellate(generators, xmin, xmax, mesh);

    // Spew the mesh statistics.
    cout << "   num mesh nodes : " << mesh.nodes.size()/2 << endl;
    cout << "   num mesh cells : " << mesh.cells.size() << endl;
    cout << "   num mesh faces : " << mesh.faces.size() << endl;
//     cout << "Node positions: " << endl;
//     for (unsigned i = 0; i != mesh.nodes.size()/2; ++i) {
//       cout << "   Node " << i << " @ (" << mesh.nodes[2*i] << " " << mesh.nodes[2*i + 1] << ")" << endl;
//     }
//     cout << "Face node sets: " << endl;
//     for (unsigned i = 0; i != nx*nx; ++i) {
//       cout << "   FACES for mesh cell " << i << " :";
//       for (unsigned j = 0; j != mesh.cells[i].size(); ++j) cout << " " << mesh.cells[i][j];
//       cout << endl;
//     }
//     for (unsigned i = 0; i != mesh.faces.size(); ++i) {
//       double xf = 0.0, yf = 0.0;
//       cout << "   NODES for mesh face " << i << " :";
//       for (unsigned j = 0; j != mesh.faces[i].size(); ++j) {
//         unsigned k = mesh.faces[i][j];
//         cout << " " << k;
//         xf += mesh.nodes[2*k];
//         yf += mesh.nodes[2*k + 1];
//       }
//       xf /= mesh.faces[i].size();
//       yf /= mesh.faces[i].size();
//       cout << " @ (" << xf << " " << yf << ")"  << endl;
//     }

    // Now do the checks.
    POLY_CHECK(mesh.nodes.size()/2 == (nx + 1)*(nx + 1));
    POLY_CHECK(mesh.cells.size() == nx*nx);
    for (unsigned i = 0; i != nx*nx; ++i) POLY_CHECK(mesh.cells[i].size() == 4);
    POLY_CHECK(mesh.faces.size() == 2*nx*(nx + 1));
  }

  cout << "PASS" << endl;
  return 0;
}
