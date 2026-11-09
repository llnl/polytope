// RotationTests
//
// A collection of rotation fields on the unit square
//
// Flow fields:
// 1. Constant Vorticity / Rigid Body Rotation
// 2. Single "Vortex in a Box" (Bell, Colella, Glas, JCP, 1989)
// 3. Taylor-Green Vortex
// 4. 16-Vortex Deformation

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "BoostTessellator.hh"
#ifdef POLYTOPE_ENABLE_TRIANGLE
#include "TriangleTessellator.hh"
#endif
#include "Communicator.hh"
#include "polytope_test_utilities.hh"

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// runTest
// -----------------------------------------------------------------------
void runTest(Tessellator<2,double>& tessellator,
             const unsigned flowType,
             const unsigned nx) {
  POLY_CHECK(flowType >= 1 and flowType <= 4);

  // Boundary size parameters
  Boundary2D boundary;
  boundary.mCenter[0] = 0.5;
  boundary.mCenter[1] = 0.5;
  boundary.mDiff = 0.4;
  boundary.setDefaultBoundary(0);
  boundary.mClipping = false;
  Generators<2> generators(boundary);
  tessellator.setPadding(boundary.mPad);
  const auto bHigh = boundary.maxBound;
  const auto bLow = boundary.minBound;
  const auto len = bHigh - bLow;
  double dx = len.x / nx;
  generators.cartesian2D(nx, nx);

  // Test name for output
  ostringstream os;
  os << "RotationTests_" << tessellator.name() << "_" << flowType;
  string testName = os.str();

  // Time stepping and point-resizing stuff
  double dt = 1.0, Tmax = 1.0, scaleFactor = 1.0;
  std::string testString = setConfig(flowType, len, dx, dt, Tmax, scaleFactor);
  std::cout << testString << std::endl;
  // Resize the generator so we don't fling them out of the boundary
  for (auto& p : generators.mPoints) {
    p = 0.5 + (p - 0.5)*scaleFactor;
  }
  // Set constant quantized bounds for the tessellator to keep a constant domain size
  std::vector<Point2<double>> bpoints = {Point2<double>(0.0, 0.0),
                                         Point2<double>(1.0, 1.0)};
  tessellator.setQBounds(bpoints);
  // The velocity field
  auto points = extractCoords<2, double>(generators.mPoints);
  vector<Point2<double>> velocityField(points.size());

  // The initial tessellation
  unsigned step = 0;
  double time = 0.0;
  Tessellation<2,double> mesh;

  // Update the point positions and generate the mesh
  vector<Point2<double>> halfTimePositions(points.size());
  while (time < Tmax) {
    mesh = tessellator.tessellate(points);
    if (step % 5 == 0) {
      outputMesh(mesh, testName, step, time);
      cout << std::round((time/Tmax)*100) << "%" << endl;
    }
    compareArea(boundary, mesh, tessellator.quantizer(), "Test " + std::to_string(flowType));
    getVelocities(points, bHigh, bLow, flowType, velocityField);
    for (unsigned i = 0; i != points.size(); ++i) {
      halfTimePositions[i] = points[i] + dt*velocityField[i];
    }
    getVelocities(halfTimePositions, bHigh, bLow, flowType, velocityField);
    for (unsigned i = 0; i != points.size(); ++i) {
      points[i] += dt*velocityField[i];
    }
    time += dt;
    ++step;
    mesh.clear();
  }
}


// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv) {
  auto& comm = Communicator::instance();
  comm.init(argc, argv);
  const int testBegin = 1;
  const int testEnd = 5;
  int nx = 20;

  if(argc >= 2) {
    nx = std::stoi(argv[1]);
  }
  POLY_CHECK2(nx > 0, "Number of points specified must be > 0.");

#ifdef POLYTOPE_ENABLE_TRIANGLE
  {
    cout << "\nTriangle Tessellator:\n" << endl;
    TriangleTessellator tessellator;
    for (unsigned flowTest = testBegin; flowTest < testEnd; ++flowTest) {
      runTest(tessellator, flowTest, nx);
    }
  }
#endif

  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator tessellator;
    for (unsigned flowTest = testBegin; flowTest < testEnd; ++flowTest) {
      runTest(tessellator, flowTest, nx);
    }
  }

  comm.finalize();
  return 0;
}
