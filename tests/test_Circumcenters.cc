// test_Circumcenters
//
// Tests for Triangle Tessellator involving degenerate Delaunay
// 1. small-area triangle has a circumcenter approaching infinity

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "polytope_test_utilities.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// test
// -----------------------------------------------------------------------
void test(Tessellator<2,double>& tessellator) {
  const int N = 40;
  double displacement = 1.0;

  string testName = "Circumcenters_" + tessellator.name();

  vector<double> PLCpoints;
  PLCpoints.push_back(0.0);  PLCpoints.push_back(0.0);
  PLCpoints.push_back(3.0);  PLCpoints.push_back(0.0);
  PLCpoints.push_back(3.0);  PLCpoints.push_back(2.0);
  PLCpoints.push_back(0.0);  PLCpoints.push_back(2.0);

  PLC<2> boundary;
  boundary.facets.resize(4, vector<int>(2));
  for (unsigned i=0; i != 4; ++i) {
     boundary.facets[i][0] = i;
     boundary.facets[i][1] = (i+1) % 4;
  }
  
  vector<double> points;
  points.push_back(1.5);  points.push_back(1.5);
  points.push_back(0.5);  points.push_back(0.5);
  points.push_back(2.5);  points.push_back(0.5);
  points.push_back(0.5);  points.push_back(1.5);
  points.push_back(1.5);  points.push_back(1.5);
  points.push_back(2.5);  points.push_back(1.5);

  double vert1[2] = {0.5, 0.5}, vert2[2] = {2.5, 0.5}, circumcenter[2];
  for (int i=0; i != N; ++i) {
     displacement *= 0.5;
     double vert3[2] = {1.5,0.5+displacement};
     geometry::computeCircumcenter2d(vert1, vert2, vert3, circumcenter);
     cout << "Step " << i << ":" << endl;
     cout << "  New generator position  = (1.5," << 0.5+displacement << ")" << endl;
     cout << "  Degenerate tri height   = " << displacement << endl;
     cout << "  Degenerate circumcenter = (" << circumcenter[0] << ","
          << circumcenter[1] << ")" << endl;
     points[1] = 0.5 + displacement;
     Tessellation<2,double> mesh;
     tessellator.tessellate(points,PLCpoints,boundary,mesh);
     //tessellator.tessellate(points,mesh);

     outputMesh(mesh, testName, points, i);

     cout << "          Number of nodes = " << mesh.nodes.size()/2 << endl;
     cout << "          Number of faces = " << mesh.faces.size() << endl << endl;
     // cout << "  Nodes:" << endl;
     // for (int j=0; j != mesh.nodes.size()/2; ++j){
     //    cout << scientific << setprecision(numeric_limits<double>::digits)
     //         <<"    (" << mesh.nodes[2*j] << "," << mesh.nodes[2*j+1] << ")" << endl;
     // }
  }
}

// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#ifdef POLYTOPE_ENABLE_MPI
   MPI_Init(&argc, &argv);
#endif

#ifdef POLYTOPE_ENABLE_TRIANGLE
  cout << "\nTriangle Tessellator:\n" << endl;
  TriangleTessellator<double> triangle;
  test(triangle);
#endif
  
#ifdef POLYTOPE_ENABLE_BOOST
  cout << "\nBoost Tessellator:\n" << endl;  
  BoostTessellator<double> boost;
  test(boost);
#endif      

  

#ifdef POLYTOPE_ENABLE_MPI
   MPI_Finalize();
#endif
   return 0;
}
