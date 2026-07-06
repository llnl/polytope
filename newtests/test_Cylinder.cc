
#include "Boundary3D.hh"
#include "polytope_test_utilities.hh"
#include "TetgenTessellator.hh"

using namespace polytope;

int main() {

  using RealType = double;
  Boundary3D boundary;
  boundary.setDefaultBoundary(3);

  std::vector<RealType> points = {0.20, 0.20, 0.20,
                                  0.80, 0.20, 0.20,
                                  0.20, 0.80, 0.20,
                                  0.20, 0.20, 0.80,
                                  0.80, 0.80, 0.80};
  Tessellation<3, double> mesh;
  TetgenTessellator tessellator;
  std::string testName = "Cylinder_tetgen";
  int test = 1;
  tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);
  outputMesh(mesh, testName, test);
  return 0;
}
