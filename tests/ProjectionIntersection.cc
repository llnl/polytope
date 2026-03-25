// test_plot
//
// Just a quick script for doing serial and parallel tests of points for

#include <fstream>
#include <vector>
#include <string>
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
void test(const Tessellator<2, double>& tessellator,
          const std::string label) {
  const int N = 40;
  for (int j = 0; j < N; ++j) {
    const double eps = 0.005 + j*0.001;
    const double x1 = 0.5 - eps;
    const double y1 = 0.5 - eps;
    const double x2 = 0.75 - 1.1*eps;
    const double y2 = 0.25 - 1.1*eps;
    const double x3 = -1.0 + 0.1401;
    const double y3 = 0.0;

    vector<double> points;
    points.push_back( 1.0);  points.push_back( 0.0);
    points.push_back( 1.0);  points.push_back(-1.0);
    points.push_back(-1.0);  points.push_back(-1.0);
    points.push_back(-1.0);  points.push_back( 1.0);
    points.push_back( 0.0);  points.push_back( 1.0);
    points.push_back( x1 );  points.push_back( y1 );
    points.push_back( x2 );  points.push_back( y2 );
    points.push_back( x3 );  points.push_back( y3 );

    ReducedPLC<2, double> boundary;
    boundary.points.push_back( 1.0);  boundary.points.push_back( 1.0);
    boundary.points.push_back(-1.0);  boundary.points.push_back( 1.0);
    boundary.points.push_back(-1.0);  boundary.points.push_back(-1.0);
    boundary.points.push_back( 1.0);  boundary.points.push_back(-1.0);
    boundary.facets.resize(4, vector<int>(2));
    for (int i = 0; i < 4; ++i) {
      boundary.facets[i][0] = i;
      boundary.facets[i][1] = (i+1)%4;
    }

    {
      Tessellation<2,double> mesh;
      tessellator.tessellate(points, mesh);
      outputMesh(mesh, "ProjectionIntersection_" + label, points, j);
    }

    {
      Tessellation<2,double> mesh;
      tessellator.tessellate(points, boundary, mesh);
      outputMesh(mesh, "ProjectionIntersection_" + label, points, N+j);     
    }
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
  {
    TriangleTessellator<double> tessellator;
    test(tessellator, "triangle");
  }
#endif
  
#ifdef POLYTOPE_ENABLE_BOOST
  {
    BoostTessellator<double> tessellator;
    test(tessellator, "boost");
  }
#endif

  cout << "PASS" << endl;

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
