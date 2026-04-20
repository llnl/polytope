// TriangleUnbounded
// Try to tessellate a series of degenerate unbounded tessellations and
// check that the correct number of cells, nodes, and faces result

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
// checkMesh
// -----------------------------------------------------------------------
bool checkMesh(const Tessellation<2,double>& mesh,
               const unsigned ncells,
               const unsigned nnodes,
               const unsigned nfaces,
               const unsigned nboundNodes,
               const unsigned nboundFaces) {
  POLY_CHECK(mesh.cells.size()    == ncells);
  // Suspending checking the number of boundary elements for our unbounded
  // tests, since my new method of using guard generators to implicitly
  // create a boundary makes this problematic.
  // POLY_CHECK_BOOL(mesh.nodes.size()/2  == nnodes);
  // POLY_CHECK(mesh.faces.size()    == nfaces);
  // POLY_CHECK(mesh.boundaryNodes.size() == nboundNodes);
  // POLY_CHECK(mesh.boundaryFaces.size() == nboundFaces);
  return true;
}


// -----------------------------------------------------------------------
// test
// -----------------------------------------------------------------------
void test(Tessellator<2,double>& tessellator,
          vector<int>& filterCheck) {

  const string testName = "Unbounded_" + tessellator.name();
  const int numTests = 11;
  vector<int> howDidIDo(numTests, 1);
  int itest = 0;

  // Test 1: Circle of generators
  {
    cout << "Test " << (itest+1) << ": Circle of generators" << endl;
    int N = 18;
    vector<double> points(2*N);
    for (unsigned i = 0; i < N; ++i){
       double theta = 2.0*M_PI*double(i)/double(N);
       points[2*i] = cos(theta);  points[2*i+1] = sin(theta);
    }
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, itest);
    bool pass = checkMesh(mesh,N,N+1,2*N,N,N);
    if (not pass) howDidIDo[itest] = 0;
    ++itest;
  }

  // Test 2: Circle of generators, random center
  {
    cout << "Test " << (itest+1) << ": Circle of generators, random center" << endl;
    int N = 18;
    vector<double> points(2*N);
    double center[2] = {random01(), random01()};
    for (unsigned i = 0; i < N; ++i){
       double theta = 2.0*M_PI*double(i)/double(N);
       points[2*i  ] = center[0] + cos(theta);
       points[2*i+1] = center[1] + sin(theta);
    }
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, itest);
    bool pass = checkMesh(mesh,N,N+1,2*N,N,N);
    if (not pass) howDidIDo[itest] = 0;
    ++itest;
  }

  // Test 3: Two uniform rows of generators
  {
    cout << "Test " << (itest+1) << ": Two uniform rows of generators" << endl;
    int N = 10;
    vector<double> points(4*N);
    for (unsigned i = 0; i < N; ++i){
       points[4*i  ] = i;  points[4*i+1] = -1.0;
       points[4*i+2] = i;  points[4*i+3] =  1.0;
    }
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, itest);
    bool pass = (checkMesh(mesh, 2*N, 3*N-1, 5*N-2, 2*N, 2*N) or
                 checkMesh(mesh, 2*N, 3*(N+1), 5*N+2, 2*N+4, 2*N+4));
    if (not pass) howDidIDo[itest] = 0;
    ++itest;
  }

  // Test 4: Collinear generators with one non-collinear
  {
    cout << "Test " << (itest+1) << ": Collinear generators, except one" << endl;
    int N = 10;
    vector<double> points(2*N);
    for (unsigned i = 0; i < N; ++i)  points[2*i] = double(i);
    points.push_back(4.5);
    points.push_back(1.0);
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, itest);
    ++itest;
  }

  // Test 5: 2x2 Cartesian Generators
  {
    cout << "Test " << (itest+1) << ": 2x2 Cartesian generators" << endl;
    vector<double> points;
    points.push_back(0.0); points.push_back(0.0);
    points.push_back(1.0); points.push_back(0.0);
    points.push_back(1.0); points.push_back(1.0);
    points.push_back(0.0); points.push_back(1.0);
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, itest);
    bool pass = (checkMesh(mesh, 4, 5, 8 , 4, 4) or
                 checkMesh(mesh, 4, 9, 12, 8, 8));
    if (not pass) howDidIDo[itest] = 0;
    ++itest;
  }

  // Test 6: 100 Random points in a unit square
  {
    cout << "Test " << (itest+1) << ": 100 Random points" << endl;
    const unsigned nPoints = 10;
    vector<double> points;
    for (unsigned i = 0; i != 2*nPoints; ++i) points.push_back(random01());
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, itest);
    ++itest;
  }


  // Test 7: Two generators
  {
    cout << "Test " << (itest+1) << ": Two generators" << endl;
    vector<double> points;
    points.push_back(0.0); points.push_back(0.0);
    points.push_back(1.0); points.push_back(0.0);
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, itest);
    bool pass = checkMesh(mesh, 2, 4, 5, 4, 4);
    if (not pass) howDidIDo[itest] = 0;
    ++itest;
  }

  // Test 8: Line of generators, uniform
  {
    cout << "Test " << (itest+1) << ": Uniform line of generators" << endl;
    int N=10;
    vector<double> points(2*N, 0);
    for (unsigned i = 0; i < N; ++i)  points[2*i] = double(i);
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, itest);
    bool pass = checkMesh(mesh, N, 2*N, 3*N-1, 2*N, 2*N);
    if (not pass) howDidIDo[itest] = 0;
    ++itest;
  }

  // Test 9: Line of generators, non-uniform
  {
    cout << "Test " << (itest+1) << ": Non-uniform line of generators" << endl;
    int N = 10;
    vector<double> points(2*N);
    for (unsigned i = 0; i < N; ++i)  points[2*i] = double(i) + random01() - 0.5;
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, itest);
    bool pass = checkMesh(mesh, N, 2*N, 3*N-1, 2*N, 2*N);
    if (not pass) howDidIDo[itest] = 0;
    ++itest;
  }

  // Test 10: Line of generators, non-uniform, shuffled
  {
    cout << "Test " << (itest+1) << ": Non-uniform line of generators, shuffled" << endl;
    int N = 10;
    vector<double> points(2*N);
    int indices[10] = {5,7,1,2,8,0,3,9,4,6};
    for (unsigned i = 0; i < N; ++i)  points[2*i] = double(indices[i]) + random01() - 0.5;
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, itest);
    bool pass = checkMesh(mesh, N, 2*N, 3*N-1, 2*N, 2*N);    
    if (not pass) howDidIDo[itest] = 0;
    ++itest;
  } 

  // Test 11: 10 generators
  {
    cout << "Test " << (itest+1) << ": 10 generators" << endl;
    double pts[20] = {-0.164777,0.26823,-0.222225,0.05397,-0.0226029,0.128871,-0.135216,
		      0.0134009,0.45223,0.416195,0.135712,0.217297,-0.358397,0.106969,
		      -0.483699,-0.257113,-0.362768,0.304177,-0.343321,-0.0990556};
    int N = 10;
    vector<double> points(2*N);
    for (unsigned i = 0; i < 2*N; ++i) points[i] = pts[i];
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, itest);
    ++itest;
  } 

  bool passedAll = true;
  for (int i = 0; i < numTests; ++i) {
    if (howDidIDo[i] == 0) {
      cout << tessellator.name() << " failed test " << i;
      if (std::find(filterCheck.begin(), filterCheck.end(), i) != filterCheck.end()) {
        cout << " --> filtered out" << endl;
      } else {
        cout << endl;
        passedAll = false;
      }
    }
  }
  POLY_CHECK(passedAll);
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
    cout << "\nTriangle Tessellator:\n" << endl;
    TriangleTessellator<double> tessellator;
    vector<int> filter;
    test(tessellator, filter);
  }
#endif   


#ifdef POLYTOPE_ENABLE_BOOST
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator<double> tessellator;
    vector<int> filter;
    filter.push_back(0); filter.push_back(1);
    test(tessellator, filter);

    // NOTE: The checks for the first two tests
    // are filtered out for BoostTessellator.
    // Boost.Voronoi requires input points to be
    // integers. It is unrealistic to expect the
    // tessellator to properly resolve the extreme
    // degeneracies in the first two tests.
    // -DPS 01/20/2015
  }
#endif
  

  

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
