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
#include "DistributedTessellator.hh"

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// computeConstantVorticityFlow
// -----------------------------------------------------------------------
void computeConstantVorticityFlow(const vector<double>& points,
                                  vector<double>& velocities) {
  const unsigned numGenerators = points.size()/2;
  for (unsigned i = 0; i != numGenerators; ++i) {
    velocities[2*i  ] = (0.5 - points[2*i+1]) * M_PI / 314;
    velocities[2*i+1] = (points[2*i  ] - 0.5) * M_PI / 314;
  }
}

// -----------------------------------------------------------------------
// computeSingleVortexFlow
// -----------------------------------------------------------------------
void computeSingleVortexFlow(const vector<double>& points,
                             vector<double>& velocities) {
  const unsigned numGenerators = points.size()/2;
  for (unsigned i = 0; i != numGenerators; ++i) {
    velocities[2*i  ] =  sin(M_PI*points[2*i]) * cos(M_PI*points[2*i+1]);
    velocities[2*i+1] = -cos(M_PI*points[2*i]) * sin(M_PI*points[2*i+1]);
  }
}

// -----------------------------------------------------------------------
// computeTaylorGreenVortexFlow
// -----------------------------------------------------------------------
void computeTaylorGreenVortexFlow(const vector<double>& points,
                                  vector<double>& velocities) {
  const unsigned numGenerators = points.size()/2;
  for (unsigned i = 0; i != numGenerators; ++i) {
    velocities[2*i  ] =  0.5 * (cos(2*M_PI*(points[2*i  ]-0.25)) *
                                sin(2*M_PI*(points[2*i+1]-0.25)) );
    velocities[2*i+1] = -0.5 * (sin(2*M_PI*(points[2*i  ]-0.25)) *
                                cos(2*M_PI*(points[2*i+1]-0.25)) );
  }
}

// -----------------------------------------------------------------------
// computeDeformationFlow
// -----------------------------------------------------------------------
void computeDeformationFlow(const vector<double>& points,
                             vector<double>& velocities) {
  const unsigned numGenerators = points.size()/2;
  double x,y;
  for (unsigned i = 0; i != numGenerators; ++i) {
    x = points[2*i  ];
    y = points[2*i+1];
    velocities[2*i  ] = -sin(4*M_PI*(x + 0.5)) * sin(4*M_PI*(y + 0.5));
    velocities[2*i+1] = -cos(4*M_PI*(x + 0.5)) * cos(4*M_PI*(y + 0.5));
  }
}

// -----------------------------------------------------------------------
// getVelocities
// -----------------------------------------------------------------------
void getVelocities(const vector<double>& points,
                   const unsigned flowType,
                   vector<double>& velocities) {
  POLY_CHECK(points.size() == velocities.size());
  switch(flowType){
  case 1:
    computeConstantVorticityFlow(points, velocities);
    break;
  case 2:
    computeSingleVortexFlow(points, velocities);
    break;
  case 3:
    computeTaylorGreenVortexFlow(points, velocities);
    break;
  case 4:
    computeDeformationFlow(points, velocities);
    break;
  }
}


// -----------------------------------------------------------------------
// runTest
// -----------------------------------------------------------------------
void runTest(Tessellator<2, double>& tessellator,
             const unsigned flowType) {
  POLY_CHECK(flowType >= 1 and flowType <= 4);

  // Boundary size parameters
  const unsigned nx = 50;

  // Figure out parallel configuration
  int rank = Communicator::getRank();
  int root = Communicator::getRoot();

  Boundary2D boundary;
  boundary.mCenter[0] = 0.5;
  boundary.mCenter[1] = 0.5;
  boundary.setDefaultBoundary(0);
  Generators<2> generators(boundary);
  generators.cartesian2D(nx, nx);
  auto finalRanks = generators.distributePointsAmongRanks();
  POLY_CONTRACT_VAR(finalRanks);
  auto& Q = Quantizer<2>::instance();
  auto bHigh = Q.m_xhi;
  auto bLow = Q.m_xlo;
  double dx = (bHigh[0] - bLow[0]) / nx;

  // Test name for output
  ostringstream os;
  os << "RotTest"
     << tessellator.name()
     << "_" << flowType;
  string testName = os.str();

  // Time stepping and point-resizing stuff
  double dt, Tmax, scaleFactor=1.0;
  double dtfactor = 10.0;
  switch(flowType){
  case 1:
    dt = 2.0;
    Tmax = 628.0/2;
    scaleFactor = 1.0/sqrt(2);
    if (rank == root) cout << "\nTest 1: Solid Rotation Flow\n" << endl;
    break;
  case 2:
    dt = sqrt(2)*dx;
    Tmax = 4.0;
    if (rank == root) cout << "\nTest 2: Single Vortex Flow\n" << endl;
    break;
  case 3:
    dt = sqrt(2)*dx;
    Tmax = 4.0;
    if (rank == root) cout << "\nTest 3: Taylor-Green (4-Vortex) Flow\n" << endl;
    break;
  case 4:
    dt = 0.5*dx;
    Tmax = 2.0;
    scaleFactor = 1.0/sqrt(2);
    if (rank == root) cout << "\nTest 4: Deformation (16-Vortex) Flow\n" << endl;
    break;
  }
  dt *= dtfactor;
  // Resize the generator so we don't fling them out of the boundary
  for (unsigned i = 0; i != generators.mPoints.size(); ++i) {
    generators.mPoints[i] = 0.5 + (generators.mPoints[i]-0.5)*scaleFactor;
  }

  // The velocity field
  vector<double> velocityField(generators.mPoints.size());

  // The initial tessellation
  unsigned step = 0;
  double time = 0.0;
  Tessellation<2,double> mesh;
  tessellator.tessellate(generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh);
  outputMesh(mesh, testName, step, time);

  // Update the point positions and generate the mesh
  vector<double> halfTimePositions(generators.mPoints.size());
  while (time < Tmax) {
    if (step % 5 == 0 and rank == root) cout << (time/Tmax)*100 << "%" << endl;
    mesh.clear();
    mesh.neighborDomains.clear();
    mesh.sharedNodes.clear();
    mesh.sharedFaces.clear();
    getVelocities(generators.mPoints, flowType, velocityField);
    for (unsigned i = 0; i != generators.mPoints.size(); ++i) {
      halfTimePositions[i] = generators.mPoints[i] + dt*velocityField[i];
    }
    getVelocities(halfTimePositions, flowType, velocityField);
    for (unsigned i = 0; i != generators.mPoints.size(); ++i) {
      generators.mPoints[i] += dt*velocityField[i];
    }
    time += dt;
    ++step;
    tessellator.tessellate(generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, testName, step, time);

    // Check the correctness of the parallel data structures
    //const string parCheck = checkDistributedTessellation(mesh);
    //POLY_CHECK2(parCheck == "ok", parCheck);
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
    // Seed the random number generator the same on all processes.
    srand(10489592);
    if (rank == root) cout << "\nTriangle Tessellator:\n" << endl;
    TriangleTessellator serialTessellator;
    DistributedTessellator<2> tessellator(serialTessellator);
    for (int i = testBegin; i < testEnd; ++i) runTest(tessellator, i);
  }
#endif

  {
    srand(10489592);
    if (rank == root) cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator serialTessellator;
    DistributedTessellator<2> tessellator(serialTessellator);
    for (int i = testBegin; i < testEnd; ++i) runTest(tessellator,i);
  }

  comm.finalize();
  return 0;
}
