#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>

#include "polytope.hh"
#include "polytope_test_utilities.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

//------------------------------------------------------------------------
// test1
//------------------------------------------------------------------------
void test1(Tessellator<2,double>& tessellator) {
  
  // Output name
  string testName = "BoostTessellator";

  // Create the generators.
  vector<double> points;
  const int nx = 3;
  const double x1 = 0.0, y1 = 0.0;
  const double x2 = 1.0, y2 = 1.0;
  const double dx = (x2 - x1)/nx, dy = (y2 - y1)/nx;
  double low [2] = {-1.0, -1.0}; // { numeric_limits<double>::max(),  numeric_limits<double>::max()};
  double high[2] = {1.0, 1.0}; // {-numeric_limits<double>::max(), -numeric_limits<double>::max()};
  unsigned ix, iy;
  double xi, yi;
  for (iy = 0; iy != nx; ++iy) {
    yi = y1 + (iy + 0.5)*dy;
    for (ix = 0; ix != nx; ++ix) {
      xi = x1 + (ix + 0.5)*dx;
      points.push_back(xi);
      points.push_back(yi);
      low [0] = min(low [0], xi);
      low [1] = min(low [1], yi);
      high[0] = max(high[0], xi);
      high[1] = max(high[1], yi);
    }
  }

  vector<double> PLCpoints;
  PLCpoints.push_back(x1);  PLCpoints.push_back(y1);
  PLCpoints.push_back(x2);  PLCpoints.push_back(y1);
  PLCpoints.push_back(x2);  PLCpoints.push_back(y2);
  PLCpoints.push_back(x1);  PLCpoints.push_back(y2);

  PLC<2,double> boundary;
  int nSides = 4;
  boundary.facets.resize(nSides, vector<int>(2));
  for (int i = 0; i != nSides; ++i) {
    boundary.facets[i][0] = i;
    boundary.facets[i][1] = (i+1)%nSides;
  }

  Tessellation<2,double> mesh;

  // Tessellate unbounded
  tessellator.tessellate(points, mesh);
  outputMesh(mesh,testName,points,0);
  mesh.clear();

  // Tessellate bounded
  tessellator.tessellate(points, PLCpoints, boundary, mesh);
  outputMesh(mesh,testName,points,1);
}


//------------------------------------------------------------------------
// test2
//------------------------------------------------------------------------
void test2(Tessellator<2,double>& tessellator) {
  
  // Output name
  string testName = "BoostTessellator";

  // Rectangular boundary
  vector<double> PLCpoints;
  PLCpoints.push_back(0.0);  PLCpoints.push_back(0.0);
  PLCpoints.push_back(1.0);  PLCpoints.push_back(0.0);
  PLCpoints.push_back(1.0);  PLCpoints.push_back(0.5);
  PLCpoints.push_back(0.0);  PLCpoints.push_back(0.5);

  // Three generators
  vector<double> points;
  points.push_back(0.2);  points.push_back(0.25   );
  points.push_back(0.6);  points.push_back(0.25   );
  points.push_back(0.8);  points.push_back(0.25001);

  // Facets
  PLC<2,double> boundary;
  int nSides = 4;
  boundary.facets.resize(nSides, vector<int>(2));
  for (int i = 0; i != nSides; ++i) {
    boundary.facets[i][0] = i;
    boundary.facets[i][1] = (i+1)%nSides;
  }
  
  Tessellation<2,double> mesh;

  // Tessellate unbounded
  tessellator.tessellate(points, mesh);
  outputMesh(mesh,testName,points,2);
  mesh.clear();

  // Tessellate bounded
  tessellator.tessellate(points, PLCpoints, boundary, mesh);
  outputMesh(mesh,testName,points,3);
}



//------------------------------------------------------------------------
// test3
//------------------------------------------------------------------------
void test3(Tessellator<2,double>& tessellator) {
  
  // Output name
  string testName = "BoostTessellator";

  // Constants
  const double degToRad = 2.0*M_PI/360.0;

  // Parameters
  const int nx = 10;
  const double degrees = 2.5;
  const double x1 = 0.0, y1 = 0.0;
  const double x2 = 1.0, y2 = 1.0;
  
  // Derived params
  const double angle = degrees * degToRad;
  const double dx = (x2 - x1)/nx, dy = (y2 - y1)/nx;
  const double xcen = 0.5*(x1+x2), ycen = 0.5*(y1+y2);
  const double cosAngle = cos(angle);
  const double sinAngle = sin(angle);

  // Create the generators.
  vector<double> points;
  double xi, yi;
  for (int iy = 0; iy != nx; ++iy) {
    yi = y1 + (iy + 0.5)*dy - ycen;
    for (int ix = 0; ix != nx; ++ix) {
      xi = x1 + (ix + 0.5)*dx - xcen;
      points.push_back(xcen + cosAngle*xi - sinAngle*yi);
      points.push_back(ycen + cosAngle*yi + sinAngle*xi);
    }
  }

  // Create the boundary
  ReducedPLC<2, double> plc;
  plc.points.push_back(x1 - 4.0*dx);  plc.points.push_back(y1 - 4.0*dy);
  plc.points.push_back(x2 + 4.0*dx);  plc.points.push_back(y1 - 4.0*dy);
  plc.points.push_back(x2 + 4.0*dx);  plc.points.push_back(y2 + 4.0*dy);
  plc.points.push_back(x1 - 4.0*dx);  plc.points.push_back(y2 + 4.0*dy);

  plc.facets.resize(4, vector<int>(2));
  for (int i = 0; i != 4; ++i) {
    plc.facets[i][0] = i;
    plc.facets[i][1] = (i+1)%4;
  }

  // The mesh
  Tessellation<2,double> mesh;

  // Tessellate unbounded
  tessellator.tessellate(points, mesh);
  outputMesh(mesh,testName,points,4);
  mesh.clear();
  
  // Tessellate bounded
  tessellator.tessellate(points, plc.points, plc, mesh);
  outputMesh(mesh,testName,points,5);
}


//------------------------------------------------------------------------
// test4
//------------------------------------------------------------------------
void test4(Tessellator<2,double>& tessellator) {
  
  // Output name
  string testName = "BoostTessellator";

  // Constants
  const double degToRad = 2.0*M_PI/360.0;

  // Parameters
  const int nx = 20;
  const double degrees = 2.5;
  const double x1 = 0.0, y1 = 0.0;
  const double x2 = 1.0, y2 = 1.0;
  
  // Derived params
  const double angle = degrees * degToRad;
  const double dx = (x2 - x1)/nx, dy = (y2 - y1)/nx;
  const double xcen = 0.5*(x1+x2), ycen = 0.5*(y1+y2);
  const double cosAngle = cos(angle);
  const double sinAngle = sin(angle);

  // PLC coordinates.
  const double xbc1 = x1 - 4.0*dx;
  const double ybc1 = y1 - 4.0*dx;
  const double xbc2 = x2 + 4.0*dx;
  const double ybc2 = y2 + 4.0*dx;
  const double xhole1 = x1 + 4.0*dx;
  const double yhole1 = y1 + 4.0*dx;
  const double xhole2 = x2 - 4.0*dx;
  const double yhole2 = y2 - 4.0*dx;
  
  // Create the generators.
  vector<double> points;
  double xi, yi, xii, yii;
  for (int iy = 0; iy != nx; ++iy) {
    yi = y1 + (iy + 0.5)*dy - ycen;
    for (int ix = 0; ix != nx; ++ix) {
      xi = x1 + (ix + 0.5)*dx - xcen;
      xii = xcen + cosAngle*xi - sinAngle*yi;
      yii = ycen + cosAngle*yi + sinAngle*xi;
      if (xii < xhole1 or xii > xhole2 or
          yii < yhole1 or yii > yhole2) {
        points.push_back(xii);
        points.push_back(yii);
      }
    }
  }

  // Create the boundary
  ReducedPLC<2, double> plc;
  plc.points.push_back(xbc1);  plc.points.push_back(ybc1);
  plc.points.push_back(xbc2);  plc.points.push_back(ybc1);
  plc.points.push_back(xbc2);  plc.points.push_back(ybc2);
  plc.points.push_back(xbc1);  plc.points.push_back(ybc2);

  plc.points.push_back(xhole1);  plc.points.push_back(yhole1);
  plc.points.push_back(xhole1);  plc.points.push_back(yhole2);
  plc.points.push_back(xhole2);  plc.points.push_back(yhole2);
  plc.points.push_back(xhole2);  plc.points.push_back(yhole1);

  plc.facets.resize(4, vector<int>(2));
  for (int i = 0; i != 4; ++i) {
    plc.facets[i][0] = i;
    plc.facets[i][1] = (i+1)%4;
  }

  plc.holes.resize(1);
  plc.holes[0].resize(4, vector<int>(2));
  for (int i = 0; i != 4; ++i) {
    plc.holes[0][i][0] = 4 + i;
    plc.holes[0][i][1] = 4 + (i+1)%4;
  }

  // The mesh
  Tessellation<2,double> mesh;

  // Tessellate unbounded
  tessellator.tessellate(points, mesh);
  outputMesh(mesh,testName,points,6);
  mesh.clear();
  
  // Tessellate bounded
  tessellator.tessellate(points, plc.points, plc, mesh);
  outputMesh(mesh,testName,points,7);
}


//------------------------------------------------------------------------
// main
//------------------------------------------------------------------------
int 
main(int argc, char** argv) 
{
#ifdef POLYTOPE_ENABLE_MPI
  MPI_Init(&argc, &argv);
#endif

  BoostTessellator<double> tessellator;
  
  {
    cout << "\nTest 1" << endl;
    test1(tessellator);
  }

  {
    cout << "\nTest 2" << endl;
    test2(tessellator);
  }

  {
    cout << "\nTest 3" << endl;
    test3(tessellator);
  }

  {
    cout << "\nTest 4" << endl;
    test4(tessellator);
  }


  

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
//------------------------------------------------------------------------
