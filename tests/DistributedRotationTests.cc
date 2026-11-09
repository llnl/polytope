// DistributedRotationTests
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
#include "polytope_test_utilities.hh"
#include "Boundary2D.hh"
#include "BoostTessellator.hh"
#ifdef POLYTOPE_ENABLE_TRIANGLE
#include "TriangleTessellator.hh"
#endif
#include "Communicator.hh"
#include "QuasiVoronoiPartitioner.hh"
#include "DistributedTessellator.hh"

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// runTest
// -----------------------------------------------------------------------
void runTest(Tessellator<2, double>& serialTessellator,
             const unsigned flowType,
             const bool cartesian = true) {
  POLY_CHECK(flowType >= 1 and flowType <= 4);

  // Boundary size parameters
  const unsigned nx = 50;

  // Figure out parallel configuration
  int rank = Communicator::getRank();
  int root = Communicator::getRoot();

  Boundary2D boundary;
  boundary.mCenter[0] = 0.5;
  boundary.mCenter[1] = 0.5;
  boundary.mDiff = 0.4;
  boundary.setDefaultBoundary(0);
  boundary.mClipping = false;
  serialTessellator.setPadding(boundary.mPad);
  Generators<2> generators(boundary);
  if (cartesian) {
    generators.cartesian2D(nx, nx);
  } else {
    const int seed = 10489592;
    const int Ngen = nx*nx;
    generators.randomPoints(Ngen, seed);
  }
  QuasiVoronoiPartitioner<2> qvpart(0);
  const auto bHigh = boundary.maxBound;
  const auto bLow = boundary.minBound;
  const auto len = bHigh - bLow;
  double dx = len.x / nx;

  // Test name for output
  ostringstream os;
  os << "RotTest"
     << serialTessellator.name()
     << "_" << flowType;
  string testName = os.str();

  // Time stepping and point-resizing stuff
  double dt = 1., Tmax = 1., scaleFactor = 1.;
  std::string testString = setConfig(flowType, len, dx, dt, Tmax, scaleFactor);
  if (rank == root) std::cout << testString << std::endl;
  // Resize the generator so we don't fling them out of the boundary
  for (auto& p : generators.mPoints) {
    p = 0.5 + (p - 0.5)*scaleFactor;
  }
  auto localPoints = qvpart.computeLocalPartition(generators.mPoints);
  // Set constant quantized bounds for the tessellator to keep a constant domain size
  std::vector<Point2<double>> bpoints = {Point2<double>(0.0, 0.0),
                                         Point2<double>(1.0, 1.0)};
  serialTessellator.setQBounds(bpoints);
  DistributedTessellator<2> tessellator(serialTessellator);

  // The velocity field
  vector<Point2<double>> velocityField(localPoints.size());

  // The initial tessellation
  unsigned step = 0;
  double time = 0.0;
  Tessellation<2,double> mesh;

  // Update the point positions and generate the mesh
  vector<Point2<double>> halfTimePositions(localPoints.size());
  while (time < Tmax) {
    mesh = tessellator.tessellate(localPoints);
    if (step % 5 == 0) {
      outputMesh(mesh, testName, step, time);
      if (rank == root) cout << std::round((time/Tmax)*100) << "%" << endl;
    }
    compareArea(boundary, mesh, tessellator.quantizer(), "Test " + std::to_string(flowType) +
                " time step " + std::to_string(step));
    getVelocities(localPoints, bHigh, bLow, flowType, velocityField);
    for (unsigned i = 0; i != localPoints.size(); ++i) {
      halfTimePositions[i] = localPoints[i] + dt*velocityField[i];
    }
    getVelocities(halfTimePositions, bHigh, bLow, flowType, velocityField);
    for (unsigned i = 0; i != localPoints.size(); ++i) {
      localPoints[i] += dt*velocityField[i];
    }
    time += dt;
    ++step;
    mesh.clear();
    Communicator::Barrier();
  }
}


// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv) {

  auto& comm = Communicator::instance();
  comm.init(argc, argv);
  int rank = Communicator::getRank();
  int root = Communicator::getRoot();

  const int testBegin = 1;
  const int testEnd   = 5;

#ifdef POLYTOPE_ENABLE_TRIANGLE
  {
    if (rank == root) cout << "\nTriangle Tessellator:\n" << endl;
    TriangleTessellator serialTessellator;
    for (int i = testBegin; i < testEnd; ++i) runTest(serialTessellator, i);
  }
#endif

  {
    if (rank == root) cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator serialTessellator;
    for (int i = testBegin; i < testEnd; ++i) runTest(serialTessellator,i);
  }

  comm.finalize();
  return 0;
}
