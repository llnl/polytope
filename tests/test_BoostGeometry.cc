#include <iostream>
#include <vector>

#include "polytope.hh"
#include "PLC_Boost_2d.hh"
#include "polytope_test_utilities.hh"

#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/io/wkt/wkt.hpp>



#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

#define BOOST_GEOMETRY_NO_ROBUSTNESS true

using namespace std;
using namespace polytope;


//------------------------------------------------------------------------------
// Convert a set of (x,y,x,y,...) coordinates to a formatted string to call
// Boost.Geometry directly
//------------------------------------------------------------------------------
template<typename IntType>
void convertToPolygonString(const IntType* pts,
			    const int nPts,
			    std::string& pstring) {
  std::ostringstream os;
  os << "POLYGON((" << pts[0] << " " << pts[1];
  for (int i = 1; i < nPts; ++i) {
    os << ", " << pts[2*i] << " " << pts[2*i+1];
  }
  os << "))";
  pstring = os.str();
}  

//------------------------------------------------------------------------------
// Intersection test - Use Polytope PLC-interface to Boost.Geometry
//------------------------------------------------------------------------------
template<typename IntType>
bool intersector_boostPLC(const IntType* pts1,
			  const int nPts1,
			  const IntType* pts2,
			  const int nPts2,
			  const bool verbose = false) {
  typedef typename polytope::ReducedPLC<2, IntType> PolygonType;
  PolygonType plc1, plc2;
  std::copy(pts1, pts1 + 2*nPts1, std::back_inserter(plc1.points));
  std::copy(pts2, pts2 + 2*nPts2, std::back_inserter(plc2.points));
  plc1.facets.resize(nPts1, std::vector<int>(2));
  plc2.facets.resize(nPts2, std::vector<int>(2));
  for (int i = 0; i < nPts1; ++i) {
    plc1.facets[i][0] = i;
    plc1.facets[i][1] = (i+1) % nPts1;
  }
  for (int i = 0; i < nPts2; ++i) {
    plc2.facets[i][0] = i;
    plc2.facets[i][1] = (i+1) % nPts2;
  }
    std::vector<PolygonType> result = polytope::BG::boost_intersect<IntType>(plc1, plc2);
  POLY_CHECK2(result.size() == 1, result.size());
  if (verbose) {
    std::cout << "Input Polygon 1: " << plc1 << std::endl;
    std::cout << "Input Polygon 2: " << plc2 << std::endl;
    std::cout << "Result: " << result[0] << std::endl;
  }
  return polytope::BG::boost_intersects<IntType>(result[0]);
}

//------------------------------------------------------------------------------
// Intersection test - Use Boost.Geometry structs directly
//------------------------------------------------------------------------------
template<typename IntType>
bool intersector_boostDirect(const std::string& s1, 
			     const std::string& s2,
			     const bool verbose = false) {
  typedef typename polytope::Point2<IntType> PointType;
  typedef boost::geometry::model::polygon<PointType, false> PolygonType;
  typedef boost::geometry::model::ring<PointType, false> RingType;
  PolygonType p1, p2;
  boost::geometry::read_wkt(s1, p1);
  boost::geometry::read_wkt(s2, p2);
  boost::geometry::correct(p1);
  boost::geometry::correct(p2);
  POLY_CHECK(not boost::geometry::intersects(p1));
  POLY_CHECK(not boost::geometry::intersects(p2));
  std::vector<PolygonType> result;
  boost::geometry::intersection(p1, p2, result);
  POLY_CHECK2(result.size() == 1, result.size());
  boost::geometry::correct(result[0]);
  if (verbose) {
    std::cout << "Input Polygon 1: " << s1 << std::endl;
    std::cout << "Input Polygon 2: " << s2 << std::endl;
    std::cout << "Result: " << std::endl;
    for (typename RingType::const_iterator itr = result[0].outer().begin();
	 itr != result[0].outer().end();
	 ++itr) std::cout << *itr << std::endl;
  }
  return boost::geometry::intersects(result[0]);
}

//------------------------------------------------------------------------------
// Driver function for the intersection tests
//------------------------------------------------------------------------------
template<typename IntType>
bool intersector(const IntType* pts1,
		 const int nPts1,
		 const IntType* pts2,
		 const int nPts2,
		 const bool useBoostPLC,
		 const bool verbose = false) {
  POLY_CHECK(nPts1 > 2 and nPts2 > 2);
  if (useBoostPLC) {
    return intersector_boostPLC<IntType>(pts1, nPts1, pts2, nPts2, verbose);
  } else {
    std::string p1str, p2str;
    convertToPolygonString<IntType>(pts1, nPts1, p1str);
    convertToPolygonString<IntType>(pts2, nPts2, p2str);
    return intersector_boostDirect<IntType>(p1str, p2str, verbose);
  }
}

template<typename IntType>
void testIntersections(const IntType* pts1,
		       const int nPts1,
		       const IntType* pts2,
		       const int nPts2,
		       const bool verbose = false) {
  std::cout << "64-bit intersection test" << std::endl;
  int64_t pts1_64[2*nPts1], pts2_64[2*nPts2];
  for (int i = 0; i < 2*nPts1; ++i) pts1_64[i] = (int64_t)(pts1[i]);
  for (int i = 0; i < 2*nPts2; ++i) pts2_64[i] = (int64_t)(pts2[i]);
  
  { // Test A - PLC interface
    std::cout << "  PLC interface:" << std::endl;
    const bool t = intersector<int64_t>(pts1_64, nPts1, pts2_64, nPts2, true, verbose);
    std::cout << "  " << (t ? "FAIL" : "PASS") << std::endl;
  }
  
  { // Test B - Boost.Geoemtry directly
    std::cout << "  Boost.Geometry directly:" << std::endl;
    const bool t = intersector<int64_t>(pts1_64, nPts1, pts2_64, nPts2, false, verbose);
    std::cout << "  " << (t ? "FAIL" : "PASS") << std::endl;
  }


  std::cout << "32-bit intersection test" << std::endl;
  int32_t pts1_32[2*nPts1], pts2_32[2*nPts2];
  for (int i = 0; i < 2*nPts1; ++i) pts1_32[i] = (int32_t)(pts1[i]);
  for (int i = 0; i < 2*nPts2; ++i) pts2_32[i] = (int32_t)(pts2[i]);  

  { // Test1a - PLC interface
    std::cout << "  PLC interface:" << std::endl;
    const bool t = intersector<int32_t>(pts1_32, nPts1, pts2_32, nPts2, true, verbose);
    std::cout << "  " << (t ? "FAIL" : "PASS") << std::endl;
  }
  
  { // Test1b - Boost.Geoemtry directly
    std::cout << "  Boost.Geometry directly:" << std::endl;
    const bool t = intersector<int32_t>(pts1_32, nPts1, pts2_32, nPts2, false, verbose);
    std::cout << "  " << (t ? "FAIL" : "PASS") << std::endl;
  }
}


//------------------------------------------------------------------------
// main
//------------------------------------------------------------------------
int main(int argc, char** argv)
{
#ifdef POLYTOPE_ENABLE_MPI
  MPI_Init(&argc, &argv);
#endif

  bool verbose;
  if (argc > 1) {
    std::cout << "Passed a variable at the command line. Turning on verbosity" << std::endl;
    verbose = true;
  } else {
    verbose = false;
  }

  { // Test 1
    std::cout << "\nTest 1" << std::endl;
    const int nPts1 = 9;
    const int nPts2 = 4;
    int32_t p1pts[18] = {42817136,-3774506,
			 43029074,-3929862,
			 31446819,18947953,
			 30772384,19615678,
			 30101303,19612322,
			 30114725,16928001,
			 33520458,6878575,
			 35332375,2413654,
			 35725796,2024148};
    int32_t p2pts[8]  = {-33386239,-33721784,
			 33721785,-33386239,
			 33386240,33721785,
			 -33721784,33386240};
    testIntersections<int32_t>(p1pts, nPts1, p2pts, nPts2, verbose);
  }


  { // Test 2
    std::cout << "\nTest 2" << std::endl;
    typedef int64_t IntType;
    const int nPts1 = 6;
    const int nPts2 = 4;
    IntType p1pts[12] = {241774944,270613768,
			 249104265,268648604,
			 258263822,277993222,
			 262822082,304885727,
			 253206525,344982011,
			 240074874,272654983};
    IntType p2pts[8]  = {225182076,225182076,
			 292290940,225182076,
			 292290940,292290940,
			 225182076,292290940};
    testIntersections<IntType>(p1pts, nPts1, p2pts, nPts2, verbose);
  }


  std::cout << "PASS" << std::endl;

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
//------------------------------------------------------------------------
