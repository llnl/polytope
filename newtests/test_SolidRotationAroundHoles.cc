// SolidRotationAroundHoles
//
// Rigid body rotation of a collection of points in a circular boundary
// with a star-shaped hole in the interior. Uses the MeshEditor to clean
// small edges as they appear as well as the cell adoption algorithm to
// handle cell pieces orphaned by boundary intersections


#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "polytope_test_utilities.hh"
#include "BoostTessellator.hh"

#include "Communicator.hh" 

#ifdef POLYTOPE_ENABLE_TRIANGLE
#include "TriangleTessellator.hh"
#endif

using namespace std;
using namespace polytope;


// -----------------------------------------------------------------------
// computeConstantVorticityFlow
// -----------------------------------------------------------------------
void computeConstantVorticityFlow(const vector<double>& points,
                                  vector<double>& velocities) {
   const unsigned numGenerators = points.size()/2;
   for (unsigned i = 0; i != numGenerators; ++i) {
      velocities[2*i  ] = -points[2*i+1] * M_PI / 314;
      velocities[2*i+1] =  points[2*i  ] * M_PI / 314;
   }
}


// -----------------------------------------------------------------------
// test
// -----------------------------------------------------------------------
void test(Tessellator<2,double>& tessellator) {
  // Initialize star-hole boundary, tessellator, and generator set
  Boundary2D boundary;
  boundary.setDefaultBoundary(5);
  vector<double> points;

  // Test name for output
  string testName = "SolidRotationAroundHoles_" + tessellator.name();

  // Timestepping parameters
  const double dt = 1.0;
  const double Tmax = 628.0/2;

  // Boundary parameters
  const double boundRadius = 1.0;
  const double theta0 = 2*M_PI/5;
  const double outerRadius = 0.75;
  const double innerRadius = outerRadius*(sin(theta0/4.0) / sin(3*theta0/4.0));

  // Add star-shaped-hole points as generators
  const bool addBoundaryGenerators = true;

  // Initialize a mask to determine if a generator will
  // move CCW (1), CW(-1), or stay fixed (0)
  vector<int> velMask;

  // Add generators between boundRadius and outerRadius
  const unsigned numRows = 6;
  const unsigned nArcs   = 90;
  const unsigned maxArc  = nArcs - 45;
  const bool alternateFlowDirection = true;
  const double dr = (boundRadius - outerRadius) / numRows;
  for (unsigned i = 0; i != numRows; ++i) {
    double r = outerRadius + (i+0.5)*dr;
    // unsigned nArcs = 6*i;
    for (unsigned j = 0; j != maxArc; ++j) {
      int direction = 1;
      if (alternateFlowDirection) direction -= 2*(i%2);
      double theta = 2*M_PI*double(j)/double(nArcs);
      double x = boundary.mCenter[0] + r*cos(theta);
      double y = boundary.mCenter[1] + r*sin(theta);
      points.push_back(x);
      points.push_back(y);
      velMask.push_back(direction);
      velMask.push_back(direction);
    }
  }
  if (addBoundaryGenerators) {
    for (unsigned j = 0; j != 5; ++j) {
      unsigned index = boundary.mPLCpoints.size()/2 - 10 + 2*j;
      points.push_back(boundary.mPLCpoints[2*index  ]);
      points.push_back(boundary.mPLCpoints[2*index+1]);
      velMask.push_back(0);
      velMask.push_back(0);
    }
  }
  // Add additional generators that will remain fixed.
  double r = 0.5*(innerRadius + outerRadius);
  unsigned numFixed = 4;
  for (unsigned i = 0; i != numFixed; ++i) {
    double theta = M_PI/2 + (i + 3.5)*theta0;
    double x = boundary.mCenter[0] + r*cos(theta);
    double y = boundary.mCenter[0] + r*sin(theta);
    points.push_back(x);
    points.push_back(y);
    velMask.push_back(0);
    velMask.push_back(0);
  }

  // The velocity field
  vector<double> velocityField(points.size());
  POLY_CHECK(velMask.size() == points.size());

  // The initial tessellation
  unsigned step = 0;
  double time = 0.0;
  Tessellation<2,double> mesh;
  tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);
  outputMesh(mesh, testName, step, time);

  // Update the point positions and generate the mesh
  vector<double> halfTimePositions(points.size());
  while (time < Tmax) {
    cout << time << endl;
    mesh.clear();
    computeConstantVorticityFlow(points, velocityField);
    for (unsigned i = 0; i != points.size(); ++i) {
      halfTimePositions[i] = points[i] + 0.5*dt*velocityField[i]*velMask[i];
    }
    computeConstantVorticityFlow(halfTimePositions, velocityField);
    for (unsigned i = 0; i != points.size(); ++i) {
      points[i] += dt*velocityField[i]*velMask[i];
    }
    time += dt;
    ++step;
    tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, testName, step, time);
  }
}


// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
  auto& comm = Communicator::instance();
  comm.init(argc, argv);

#ifdef POLYTOPE_ENABLE_BOOST
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator tessellator;
    test(tessellator);
  }
#endif

#ifdef POLYTOPE_ENABLE_TRIANGLE
  {
    cout << "\nTriangle Tessellator:\n" << endl;
    TriangleTessellator tessellator;
    test(tessellator);
  }
#endif

  comm.finalize();
  return 0;
}
