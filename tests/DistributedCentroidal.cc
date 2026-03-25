// test_DistributedCentroidal
//
// Perform repeated Lloyd iterations on the points inside a specified boundary
// to converge (slowly) to a centroidal Voronoi. Points are distributed across
// domain boundaries, and communication is updated with each iteration

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <algorithm>

#include "polytope.hh"
#include "polytope_test_utilities.hh"
#include "polytope_geometric_utilities.hh"

#include "mpi.h"

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
// lloyd
// -----------------------------------------------------------------------
void lloyd(Tessellation<2,double>& mesh,
           vector<double>& points) {
   for (unsigned i = 0; i < mesh.cells.size(); ++i){
      double cent[2], area;
      geometry::computeCellCentroidAndSignedArea(mesh, i, 1.0e-12, cent, area);
      points[2*i  ] = 0.5*(points[2*i  ] + cent[0]);
      points[2*i+1] = 0.5*(points[2*i+1] + cent[1]);
   }
}

// -----------------------------------------------------------------------
// conditionNumberNodes
// -----------------------------------------------------------------------
vector<double> conditionNumberNodes(Tessellation<2,double>& mesh) {
   vector<double> cond(mesh.nodes.size()/2, 0.0);
   for (unsigned icell = 0; icell != mesh.cells.size(); ++icell) {
      for (unsigned f = 0; f != mesh.cells[icell].size(); ++f) {
         int f1 = mesh.cells[icell][f];
         int f2 = mesh.cells[icell][(f+1) % mesh.cells[icell].size()];
         const unsigned iface1 = (f1 < 0) ? ~f1 : f1;
         const unsigned iface2 = (f2 < 0) ? ~f2 : f2;
         POLY_ASSERT(iface1 < mesh.faces.size() and iface2 < mesh.faces.size());
         POLY_ASSERT(mesh.faces[iface1].size() == 2 and mesh.faces[iface2].size() == 2);
         POLY_ASSERT(iface1 != iface2);
         const unsigned inode1 = (f1 >= 0) ? mesh.faces[iface1][0] : mesh.faces[iface1][1];
         const unsigned inode2 = (f1 >= 0) ? mesh.faces[iface1][1] : mesh.faces[iface1][0];
         const unsigned itmp1  = (f2 >= 0) ? mesh.faces[iface2][0] : mesh.faces[iface2][1];
         const unsigned itmp2  = (f2 >= 0) ? mesh.faces[iface2][1] : mesh.faces[iface2][0];
         const unsigned inode3 = (itmp1 == inode1 or itmp1 == inode2) ? itmp2 : itmp1;
         POLY_ASSERT(inode2 < mesh.nodes.size()/2);
         POLY_ASSERT(inode1 != inode2);
         POLY_ASSERT(inode2 != inode3);
         POLY_ASSERT(inode1 != inode3);
         const double lp = geometry::distance<2,double>(&mesh.nodes[2*inode1], &mesh.nodes[2*inode2]);
         const double lq = geometry::distance<2,double>(&mesh.nodes[2*inode2], &mesh.nodes[2*inode3]);
         const double lr = geometry::distance<2,double>(&mesh.nodes[2*inode3], &mesh.nodes[2*inode1]);
         const double P = 0.5*(lp + lq + lr);
         const double A = sqrt(P*(P-lp)*(P-lq)*(P-lr));
         cond[inode2] += (lp*lp + lq*lq) / A;
      }
   }
   return cond;
}

// -----------------------------------------------------------------------
// conditionNumber
// -----------------------------------------------------------------------
vector<double> conditionNumber(Tessellation<2,double>& mesh) {
  vector<double> cond(mesh.cells.size(), 0.0);
  for (unsigned icell = 0; icell != mesh.cells.size(); ++icell) {
    double cent[2], area;
    geometry::computeCellCentroidAndSignedArea(mesh, icell, 1.0e-12, cent, area);
    for (vector<int>::const_iterator itr = mesh.cells[icell].begin();
         itr != mesh.cells[icell].end(); ++itr) {
      const unsigned iface = (*itr < 0) ? ~(*itr) : *itr;
      POLY_ASSERT(iface < mesh.faces.size());
      POLY_ASSERT(mesh.faces[iface].size() == 2);
      const unsigned inode1 = (*itr < 0) ? mesh.faces[iface][1] : mesh.faces[iface][0];
      const unsigned inode2 = (*itr < 0) ? mesh.faces[iface][0] : mesh.faces[iface][1];
      POLY_ASSERT(inode1 != inode2);
      POLY_ASSERT(inode1 < mesh.nodes.size()/2 and inode2 < mesh.nodes.size()/2);
      const double lp = geometry::distance<2,double>(cent, &mesh.nodes[2*inode1]);
      const double lq = geometry::distance<2,double>(cent, &mesh.nodes[2*inode2]);
      const double lr = geometry::distance<2,double>(&mesh.nodes[2*inode1], &mesh.nodes[2*inode2]);
      const double P = 0.5*(lp + lq + lr);
      const double A = sqrt(P*(P-lp)*(P-lq)*(P-lr));
      cond[icell] += (lp*lp + lq*lq) / A;
    }
  }
  return cond;
}

// -----------------------------------------------------------------------
// lloydTestDistributed
// -----------------------------------------------------------------------
void lloydTestDistributed(Tessellator<2,double>& tessellator) {
  const unsigned nPoints     = 1000;     // Number of generators
  const unsigned nIter       = 2000;     // Number of iterations
  unsigned outputEvery = 10;      // Output frequency
  const int btype = 2;

  // Seed the random number generator the same on all processes.
  srand(10489592);
  
  // Test name
  string testName = "DistributedLloydTest_" + tessellator.name();

  // Figure out our parallel configuration.
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // Set up boundary and disperse random generator locations
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(btype);
  Generators<2,double> generators(boundary);
  generators.randomPoints(nPoints);

  // Assign points to processors in quasi-Voronoi fashion
  vector<double> points, xproc, yproc;
  double p[2];
  xproc.reserve(numProcs);
  yproc.reserve(numProcs);
  for (unsigned iproc = 0; iproc != numProcs; ++iproc) {
    boundary.getPointInside(p);
    xproc.push_back(p[0]);
    yproc.push_back(p[1]);
  }
  for (unsigned i = 0; i < nPoints; ++i){
    unsigned owner = 0;
    double minDist2 = distance2(generators.mPoints[2*i], 
                                generators.mPoints[2*i+1], 
                                xproc[0], yproc[0]);
    for (unsigned iproc = 1; iproc < numProcs; ++iproc) {
      const double d2 = distance2(generators.mPoints[2*i], 
                                  generators.mPoints[2*i+1], 
                                  xproc[iproc], yproc[iproc]);
      if( d2 < minDist2 ){
        owner = iproc;
        minDist2 = d2;
      }
    }
    if (rank == owner) {
      points.push_back(generators.mPoints[2*i  ]);
      points.push_back(generators.mPoints[2*i+1]);
    }
  }
  
  // Initialize mesh and tessellator
  Tessellation<2,double> mesh;
  tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);

  // Do the Lloyd iteration thang
  unsigned iter = 0;
  vector<double> cond = conditionNumber(mesh);
  outputMesh(mesh, testName, points, cond, iter);
  while (iter != nIter) {
    lloyd(mesh, points);
    ++iter;
    mesh.clear();
    tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);
    if (iter % outputEvery == 0 or iter <= 20) {
       if (rank == 0) cerr << iter << endl;
       cond = conditionNumber(mesh);
       outputMesh(mesh, testName, points, cond, iter);
    }
    if (iter > 200) outputEvery = 100;
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
   
#ifdef POLYTOPE_ENABLE_TRIANGLE
  {
    cout << "\nTriangle Tessellator:\n" << endl;
    DistributedTessellator<2,double> tessellator(new TriangleTessellator<double>(), true, true);
    lloydTestDistributed(tessellator);
  }
#endif   

#ifdef POLYTOPE_ENABLE_BOOST
  {
    cout << "\nBoost Tessellator:\n" << endl;
    DistributedTessellator<2,double> tessellator(new BoostTessellator<double>(), true, true);
    lloydTestDistributed(tessellator);
  }
#endif


  
   
  MPI_Finalize();
  return 0;
}
