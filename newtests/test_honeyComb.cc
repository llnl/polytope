// test_honeyComb

#include "polytope.hh"
#include "polytope_test_utilities.hh"
#include "BoostTessellator.hh"
#include "Generators.hh"
#include "SiloReader.hh"
#include "SiloUtils.hh"
#include "Communicator.hh"

#include <vector>

#ifdef POLYTOPE_ENABLE_TRIANGLE
#include "TriangleTessellator.hh"
#endif

using namespace polytope;

// -----------------------------------------------------------------------
// test
// -----------------------------------------------------------------------
void test(Tessellator<2,double>& tessellator, const std::string& outname, bool honeyComb) {
  double xmin = -0.5;
  double ymin = -0.5;
  double xmax = 0.5;
  double ymax = 0.5;
  unsigned nx = 8;
  unsigned ny = 8;
  std::vector<double> points;
  int cycle = 0;
  int numNodes = 0;
  if (honeyComb) {
    double dx = (xmax - xmin)/(nx + 1.);
    double dy = (ymax - ymin)/(ny + 1.);
    cycle = 1;
    int iy = -1;
    for (int ix = 0; ix <= nx; ++ix) {
      if ((ix % 2) > 0) {
        double x = xmin + (ix + 0.5) * dx;
        double y = ymin + (iy + 0.5 + 0.5*(ix % 2)) * dy;
        points.push_back(x);
        points.push_back(y);
      }
    }
    for (int iy = 0; iy <= ny; ++iy) {
      for (int ix = 0; ix <= nx; ++ix) {
        double x = xmin + (ix + 0.5) * dx;
        double y = ymin + (iy + 0.5 + 0.5*(ix % 2)) * dy;
        points.push_back(x);
        points.push_back(y);
      }
    }
    numNodes = 4 + 2*(nx + ny + nx*ny) + nx;
  } else {
    numNodes = (nx+1)*(ny+1);
    double dx = (xmax - xmin)/double(nx);
    double dy = (ymax - ymin)/double(ny);
    // Cartesian case: Just put the generators at the proposed zone centers.
    // But note that none of the generators are actually on the boundary, unlike
    // in the honeycomb case.
    for (int iy = 0; iy < ny; ++iy) {
      double y = ymin + (iy + 0.5) * dy;
      for (int ix = 0; ix < nx; ++ix) {
        double x = xmin + (ix + 0.5) * dx;
        points.push_back(x);
        points.push_back(y);
      }
    }
  }
  Point2<double> lo(xmin, ymin);
  Point2<double> hi(xmax, ymax);
  std::vector<double> plcPoints = flattenCoords(shapes::createSquarePoints(lo, hi));
  PLC<2> plc;
  plc.facets = shapes::createSquareFaces();
  auto& Q = Quantizer<2>::instance();
  Q.init(plcPoints);
  Tessellation<2, double> mesh;
  tessellator.tessellate(points, mesh);
  outputMesh(mesh, outname, cycle);
  testWatertight(mesh, 0);
  // Read the mesh we just wrote back in
  Tessellation<2, double> readMesh;
  std::string masterFilename = getMasterFilename(outname, cycle);
  SiloReader<2, Tessellation<2, double>>::FieldTypeMap fields;
  SiloReader<2, Tessellation<2, double>>::read(readMesh, fields, masterFilename);
  // Try writing the mesh back out
  outputMesh(readMesh, "re"+outname, cycle);
  int meshnodes = readMesh.nodes.size();
  POLY_CHECK2(meshnodes == numNodes, "Number of nodes in mesh: " << meshnodes
              << " does not match expected number of nodes: " << numNodes);
}
// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv) {
  auto& comm = Communicator::instance();
  comm.init(argc, argv);

#ifdef POLYTOPE_ENABLE_BOOST
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator tessellator;
    test(tessellator, "boosthoney", false);
    test(tessellator, "boosthoney", true);
  }
#endif

#ifdef POLYTOPE_ENABLE_TRIANGLE
  {
    cout << "\nTriangle Tessellator:\n" << endl;
    TriangleTessellator tessellator;
    test(tessellator, "trianglehoney", false);
    test(tessellator, "trianglehoney", true);
  }
#endif

  comm.finalize();
  return 0;
}
