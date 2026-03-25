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
#include "checkDistributedTessellation.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;


//------------------------------------------------------------------------------
// Compute the square of the distance.
//------------------------------------------------------------------------------
double distance2(const double x1, const double y1,
                 const double x2, const double y2) {
  return (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1);
}

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
   POLY_ASSERT(points.size() == velocities.size());
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
void runTest(Tessellator<2,double>& tessellator,
             const unsigned flowType) {
  POLY_ASSERT(flowType >= 1 and flowType <= 4);
  
  // Boundary size parameters
  const unsigned nx = 50;
  const double xmin = 0.0, xmax = 1.0;
  const double ymin = 0.0, ymax = 1.0;
  const double dx = (xmax-xmin)/nx,  dy = (ymax-ymin)/nx;

  // Figure out parallel configuration
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // Seed the random number generator the same on all processes.
  srand(10489592);

  // Create the seed positions for each domain.  Note we rely on this sequence
  // being the same for all processors and therefore don't need to communicate
  // this information.
  vector<double> xproc, yproc;
  xproc.reserve(numProcs);
  yproc.reserve(numProcs);
  for (unsigned iproc = 0; iproc != numProcs; ++iproc) {
    xproc.push_back(xmin + random01()*(xmax - xmin));
    yproc.push_back(ymin + random01()*(ymax - ymin));
  }

  // Test name for output
  ostringstream os;
  os << "DistributedRotationTests_" 
     << tessellator.name() 
     << "_" << flowType
     << "_" << numProcs << "procs";
  string testName = os.str();
   
  // Time stepping and point-resizing stuff
  double dt, Tmax, scaleFactor=1.0;
  switch(flowType){
  case 1:
    dt = 2.0;
    Tmax = 628.0/2;
    scaleFactor = 1.0/sqrt(2);
    if (rank == 0) cout << "\nTest 1: Solid Rotation Flow\n" << endl;
    break;
  case 2:
    dt = sqrt(2)*dx;
    Tmax = 4.0;
    if (rank == 0) cout << "\nTest 2: Single Vortex Flow\n" << endl;
    break;
  case 3:
    dt = sqrt(2)*dx;
    Tmax = 4.0;
    if (rank == 0) cout << "\nTest 3: Taylor-Green (4-Vortex) Flow\n" << endl;
    break;
  case 4:
    dt = 0.5*dx;
    Tmax = 2.0;
    scaleFactor = 1.0/sqrt(2);
    if (rank == 0) cout << "\nTest 4: Deformation (16-Vortex) Flow\n" << endl;
    break;
  }

  // The PLC points for a unit square
  vector<double> PLCpoints;
  PLCpoints.push_back(xmin);  PLCpoints.push_back(ymin);
  PLCpoints.push_back(xmax);  PLCpoints.push_back(ymin);
  PLCpoints.push_back(xmax);  PLCpoints.push_back(ymax);
  PLCpoints.push_back(xmin);  PLCpoints.push_back(ymax);

  // The unit square PLC facets
  PLC<2,double> boundary;
  boundary.facets.resize(4, vector<int>(2));
  for (int i = 0; i != 4; ++i) {
    boundary.facets[i][0] = i;
    boundary.facets[i][1] = (i+1)%4;
  }
   
  // The generator set
  vector<double> points;
  unsigned ix, iy;
  double xi, yi;
  for (iy = 0; iy != nx; ++iy) {
    yi = ymin + (iy + 0.5)*dy;
    for (ix = 0; ix != nx; ++ix) {
      xi = xmin + (ix + 0.5)*dx;
      unsigned owner = 0;
      double minDist2 = distance2(xi, yi, xproc[0], yproc[0]);
      for (unsigned iproc = 1; iproc < numProcs; ++iproc) {
        const double d2 = distance2(xi, yi, xproc[iproc], yproc[iproc]);
        if (d2 < minDist2) {
          owner = iproc;
          minDist2 = d2;
        }
      }
      if (rank == owner) {
        points.push_back(xi);
        points.push_back(yi);
      }
    }
  }

  POLY_CHECK2(points.size()/2 > 1, "Processor " << rank << " only has "
              << points.size()/2 << " generators");

  // Resize the generator so we don't fling them out of the boundary
  for (unsigned i = 0; i != points.size(); ++i) {
    points[i] = 0.5 + (points[i]-0.5)*scaleFactor;  
  }

  // The velocity field
  vector<double> velocityField(points.size());

  // The initial tessellation
  unsigned step = 0;
  double time = 0.0;
  Tessellation<2,double> mesh;
  tessellator.tessellate(points, PLCpoints, boundary, mesh);
  outputMesh(mesh, testName, points, step, time);
   
  // Update the point positions and generate the mesh
  vector<double> halfTimePositions(points.size());
  while (time < Tmax) {
    if (step % 5 == 0 and rank == 0) cout << (time/Tmax)*100 << "%" << endl;
    //if (rank == 0) cout << step << endl;
    mesh.clear();
    mesh.neighborDomains.clear();
    mesh.sharedNodes.clear();
    mesh.sharedFaces.clear();
    getVelocities(points, flowType, velocityField);
    for (unsigned i = 0; i != points.size(); ++i) {
      halfTimePositions[i] = points[i] + dt*velocityField[i];
    }
    getVelocities(halfTimePositions, flowType, velocityField);
    for (unsigned i = 0; i != points.size(); ++i) {
      points[i] += dt*velocityField[i];
      POLY_ASSERT(xmin <= points[i] and points[i] <= xmax);
    }
    time += dt;
    ++step;
    tessellator.tessellate(points, PLCpoints, boundary, mesh);
    outputMesh(mesh, testName, points, step, time);

    // Check the correctness of the parallel data structures
    const string parCheck = checkDistributedTessellation(mesh);
    POLY_ASSERT2(parCheck == "ok", parCheck);
    MPI_Barrier(MPI_COMM_WORLD);
  }
}


// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  const int testBegin = 1;
  const int testEnd   = 5;

#ifdef POLYTOPE_ENABLE_TRIANGLE
  {
    if (rank==0) cout << "\nTriangle Tessellator:\n" << endl;
    DistributedTessellator<2, double> tessellator(new TriangleTessellator<double>(),
                                                  true, true);
    for (int i = testBegin; i < testEnd; ++i) runTest(tessellator,i);
  }
#endif   


#ifdef POLYTOPE_ENABLE_BOOST
  {
    if (rank==0) cout << "\nBoost Tessellator:\n" << endl;
    DistributedTessellator<2, double> tessellator(new BoostTessellator<double>(),
                                                  true, true);
    for (int i = testBegin; i < testEnd; ++i) runTest(tessellator,i);
  }
#endif
   

  cout << "PASS" << endl;

  MPI_Finalize();
  return 0;
}
