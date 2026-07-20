// test_OrphanCases
//
// A collection of difficult test cases for the orphaned cell algorithm.
// All tests involve meshing a circular region with a star-shaped hole
// in the middle with only 20 points. Seeding the random number generator
// provides the input generator locations. Test cases were found through 
// trial-and-error (though with obnoxious frequency!).
// -----------------------------------------------------------------------

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "polytope_test_utilities.hh"
#include "GeomUtils.hh"
#ifdef POLYTOPE_ENABLE_TRIANGLE
#include "TriangleTessellator.hh"
#endif
#include "BoostTessellator.hh"

#include "Communicator.hh" 

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// printArea
// -----------------------------------------------------------------------
void printArea(Boundary2D<double>& boundary,
	       Tessellation<2,double>& mesh) {
   const double area = computeTessellationArea(mesh);
   const double relErr = std::abs(boundary.mArea-area)/boundary.mArea;
   cout << "Tessellation Area = " << area << endl;
   cout << "Relative error    = " << relErr << endl;
   POLY_CHECK(relErr < 1.0E-8);
}

// -----------------------------------------------------------------------
// checkNearestNode
// -----------------------------------------------------------------------
bool checkNearestNode(const Tessellation<2,double>& mesh,
		      const double /*tol*/) {
  set<unsigned> boundaryNodes;
  for (unsigned iface = 0; iface != mesh.faces.size(); ++iface) {
    if (mesh.faceCells[iface].size() == 1) {
      for (vector<unsigned>::const_iterator itr = mesh.faces[iface].begin();
	   itr != mesh.faces[iface].end();
	   ++itr ) {
        boundaryNodes.insert(*itr);
      }
    }
  }

  vector<double> minDistList(mesh.nodes.size()/2);
  double dist, minDist, mostMin = numeric_limits<double>::max();
  for (unsigned i = 0; i != mesh.nodes.size()/2; ++i) {
    minDist = numeric_limits<double>::max();
    for (unsigned j = 0; j != mesh.nodes.size()/2; ++j) {
      if (i != j) {
	dist = distance<2, double>(&mesh.nodes[2*i], &mesh.nodes[2*j]);
	minDist = std::min(dist, minDist);
      }
    }
    mostMin = std::min(mostMin, minDist);
    minDistList[i] = minDist;
  }

  // cerr << "Boundary node distances:" << endl;
  // for (set<unsigned>::const_iterator itr = boundaryNodes.begin();
  //      itr != boundaryNodes.end();
  //      ++itr ) cerr << "   " << "(" << mesh.nodes[2*(*itr)]
  //       	    << "," << mesh.nodes[2*(*itr)+1] << "): "
  //       	    << minDistList[*itr] << endl;

  // cerr << "Minimum node-node distance = " << mostMin << endl;
  // return (mostMin > tol) ? true : false;
  return true;
}


// -----------------------------------------------------------------------
// test
// -----------------------------------------------------------------------
void test(Tessellator<2,double>& tessellator) {

  // Test name for output
  string testName = "OrphanCases_" + tessellator.name();
  
  // Initialize boundary and tessellator
  Boundary2D<double> boundary;
  
  // Circular region with star-shaped hole
  int bType = 5;
  boundary.setDefaultBoundary(bType);
  
  int i = 1;
  const double dist = 1.0e-6;
  auto& Q = Quantizer<2>::instance();
  
  // Test 1: Cell parents multiple orphans
  int seed = 10489593;
  {
    cout << "\nTest 1: Cell parents multiple orphans" << endl;
    Generators<2,double> generators(boundary);
    generators.randomPoints(20, seed);
    Tessellation<2,double> mesh;
    tessellator.tessellate(generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, testName,i);
    printArea(boundary,mesh);
    POLY_CHECK(checkNearestNode(mesh, dist));
    ++i;
  }
  
  // Test 2: Orphan neighbors are also parents of orphans
  {
    seed++;
    cout << "\nTest 2: Orphan neighbors are also parents of orphans" << endl;
    Generators<2,double> generators(boundary);
    generators.randomPoints(20, seed);
    Tessellation<2,double> mesh;
    tessellator.tessellate(generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, testName,i);
    printArea(boundary,mesh);
    POLY_CHECK(checkNearestNode(mesh,dist));
    ++i;
  }
  
  // Test 3: Overlapping orphans
  {
    seed = 10489609;
    cout << "\nTest 3: Overlapping orphans" << endl;
    Generators<2,double> generators(boundary);
    generators.randomPoints(20, seed);
    Tessellation<2,double> mesh;
    tessellator.tessellate(generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, testName, i);
    printArea(boundary,mesh);
    POLY_CHECK(checkNearestNode(mesh, dist));
    ++i;
  }
  
  // Test 4: Empty orphan neighbor set
  {
    seed = 10489611;
    cout << "\nTest 4: Empty orphan neighbor set" << endl;
    Generators<2,double> generators(boundary);
    generators.randomPoints(20, seed);
    Tessellation<2,double> mesh;
    tessellator.tessellate(generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, testName, i);
    printArea(boundary,mesh);
    POLY_CHECK(checkNearestNode(mesh, dist));
    ++i;
  }
  
  // Test 5: Boost.Geometry calls invalid overlay exception
  {
    seed = 10489612;
    cout << "\nTest 5: Boost.Geometry calls invalid overlay exception" << endl;
    Generators<2,double> generators(boundary);
    generators.randomPoints(20, seed);
    Tessellation<2,double> mesh;
    tessellator.tessellate(generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, testName, i);
    printArea(boundary,mesh);
    POLY_CHECK(checkNearestNode(mesh, dist));
    ++i;
  }
  
  // Test 6: Nonconvex boundary with three internal generators
  {
    cout << "\nTest 6: Nonconvex boundary with three internal generators" << endl;
    std::vector<double> points;
    points.push_back(0.00); points.push_back(0.60);
    points.push_back(0.40); points.push_back(0.10);
    points.push_back(0.60); points.push_back(0.50);
    std::vector<double> PLCpoints;
    PLCpoints.push_back(0.0); PLCpoints.push_back(0.0);
    PLCpoints.push_back(0.1); PLCpoints.push_back(0.0);
    PLCpoints.push_back(0.2); PLCpoints.push_back(0.8);
    PLCpoints.push_back(0.3); PLCpoints.push_back(0.0);
    PLCpoints.push_back(1.0); PLCpoints.push_back(0.0);
    PLCpoints.push_back(1.0); PLCpoints.push_back(1.0);
    PLCpoints.push_back(0.0); PLCpoints.push_back(1.0);

    for (unsigned k = 0; k != PLCpoints.size(); ++k) {
      points.push_back(PLCpoints[k]);
    }

    PLC<2> boundary;
    boundary.facets.resize(7, std::vector<int>(2));
    for (int j = 0; j < 7; ++j){
      boundary.facets[j][0] = j;
      boundary.facets[j][1] = (j+1) % 7;
    }
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, PLCpoints, boundary, mesh);
    outputMesh(mesh, testName, i);
    ++i;
  }
  
  // Test 7: Original 3x3 Test Case
  {
    cout << "\nTest 7: 3x3 Unit Test with 2 Orphans" << endl;
    vector<double> PLCpoints;
    vector<double> points;
    PLC<2> boundary;
    Tessellation<2,double> mesh;
    PLCpoints.push_back(0.0);  PLCpoints.push_back(0.0);
    PLCpoints.push_back(1.2);  PLCpoints.push_back(0.0);
    PLCpoints.push_back(1.2);  PLCpoints.push_back(1.3);
    PLCpoints.push_back(1.3);  PLCpoints.push_back(1.3);
    PLCpoints.push_back(1.3);  PLCpoints.push_back(0.0);
    PLCpoints.push_back(3.0);  PLCpoints.push_back(0.0);
    PLCpoints.push_back(3.0);  PLCpoints.push_back(3.0);
    PLCpoints.push_back(1.3);  PLCpoints.push_back(3.0);
    PLCpoints.push_back(1.3);  PLCpoints.push_back(1.7);
    PLCpoints.push_back(1.2);  PLCpoints.push_back(1.7);
    PLCpoints.push_back(1.2);  PLCpoints.push_back(3.0);
    PLCpoints.push_back(0.0);  PLCpoints.push_back(3.0);
    
    int ix, iy, nx = 3;
    double xi, yi;
    for (iy = 0; iy != nx; ++iy) {
      yi = iy + 0.5;
      for (ix = 0; ix != nx; ++ix) {
	xi = ix + 0.5;
	points.push_back(xi);  points.push_back(yi);
      }
    }
    
    int nSides = PLCpoints.size()/2;
    boundary.facets.resize( nSides, std::vector<int>(2) );
    for (unsigned j = 0; j != nSides; ++j){
      boundary.facets[j][0] = j;
      boundary.facets[j][1] = (j+1) % nSides;
    }
    Q.init(PLCpoints, -1, 0.5);
    tessellator.tessellate(points, PLCpoints, boundary, mesh);
    outputMesh(mesh, testName, i);
    const double trueArea = 8.74;
    const double tessArea = computeTessellationArea(mesh);
    const double fracerr  = std::abs(trueArea - tessArea)/trueArea;
    POLY_CHECK2(fracerr < 1.0E-8, "Relative error in the tessellation "
                 << "area exceeds tolerance:" << endl
                 << "      Area = " << tessArea << endl
                 << "     Error = " << trueArea - tessArea << endl
                 << "Frac Error = " << fracerr);
    POLY_CHECK(checkNearestNode(mesh, dist));
    ++i;
  }

  // Test 8: Lots of random points
  {
    cout << "\nTest 8: Lots of random points" << endl;
    const unsigned N = 100;
    seed = 10332520;
    Q.init(boundary.mPLCpoints, -1., 0.1);
    for (unsigned iter = 0; iter < N; ++iter) {
      Generators<2,double> generators(boundary);
      generators.randomPoints(50, seed);
      Tessellation<2,double> mesh;
      tessellator.tessellate(generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh);
      outputMesh(mesh, testName, i+iter);
      cout << iter << endl;
      printArea(boundary,mesh);
      bool result = checkNearestNode(mesh, dist);
      POLY_CHECK(result);
      seed++;
    }
    ++i;
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
    cout << "\nTriangle Tessellator:\n" << endl;
    TriangleTessellator tessellator;
    test(tessellator);
  }
#endif   

#ifdef POLYTOPE_ENABLE_BOOST
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator tessellator;
    test(tessellator);
  }
#endif   

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Finalize();
#endif
   return 0;
}
