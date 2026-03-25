// Unit tests for our hashing algorithms.

#include <iostream>
#include <vector>
#include <map>
#include <cassert>
#include <ctime>

#include "polytope.hh"
#include "polytope_internal.hh"
#include "polytope_geometric_utilities.hh"
#include "polytope_test_utilities.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

namespace {

//------------------------------------------------------------------------------
// The common checking code for our various hashing tests.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void hashSlingingSlasher(const RealType* xmin,
                         const RealType* xmax,
                         const RealType tol,
                         const std::vector<RealType>& coords) {
  POLY_CHECK(coords.size() % Dimension == 0);
  typedef polytope::geometry::Hasher<Dimension, RealType> HasherType;
  const unsigned n = coords.size() / Dimension;
  std::map<uint64_t, std::vector<unsigned> > hash2ids;
  RealType pos[Dimension];
  for (unsigned i = 0; i != n; ++i) {
    for (unsigned j = 0; j != Dimension; ++j) {
      POLY_CHECK(coords[Dimension*i + j] >= xmin[j] and coords[Dimension*i + j] <= xmax[j]);
    }
    const uint64_t hashi = HasherType::hashPosition(&coords[Dimension*i], xmin, xmax, xmin, xmax, tol);
    hash2ids[hashi].push_back(i);
    for (std::vector<unsigned>::const_iterator itr = hash2ids[hashi].begin();
         itr != hash2ids[hashi].end();
         ++itr) {
      POLY_CHECK2((geometry::distance<Dimension, RealType>(&coords[Dimension*i], &coords[Dimension*(*itr)]) <= 2*tol),
                  i << " and " << (*itr) << " hashed to same position, but not!  "
                  << (geometry::distance<Dimension, RealType>(&coords[Dimension*i], &coords[Dimension*(*itr)])));
    }
    HasherType::unhashPosition(pos, xmin, xmax, xmin, xmax, hashi, tol);
    POLY_CHECK2((geometry::distance<Dimension, RealType>(&coords[Dimension*i], pos) <= tol),
                "Unhash position too far : " << (geometry::distance<Dimension, RealType>(&coords[Dimension*i], pos)));
    // cerr << i << " (" << coords[2*i] << " " << coords[2*i+1] << ") --> " << hashi << " : (" 
    //      << pos[0] << " " << pos[1] << ") --> " << HasherType::hashPosition(pos, xmin, xmax, xmin, xmax, tol) << endl;
    POLY_CHECK(HasherType::hashPosition(pos, xmin, xmax, xmin, xmax, tol) == hashi);
  }
}

}

// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv) {

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Init(&argc, &argv);
#endif

  const unsigned n = 1000;

  //----------------------------------------------------------------------
  // Hash a random set of points in 2D using a fine tolerance.
  //----------------------------------------------------------------------
  {
    cout << "2D random hash test, doubles, fine tolerance." << endl;
    const double xmin[2] = {-10.0, -10.0}, xmax[2] = {10.0, 10.0};
    const double tol = 1.5e-5;
    std::vector<double> coords(2*n);
    for (unsigned i = 0; i != n; ++i) {
      coords[2*i  ] = xmin[0] + random01()*(xmax[0] - xmin[0]);
      coords[2*i+1] = xmin[1] + random01()*(xmax[1] - xmin[1]);
    }
    hashSlingingSlasher<2, double>(xmin, xmax, tol, coords);
  }

  //----------------------------------------------------------------------
  // Hash a random set of points in 2D using a coarse tolerance.
  //----------------------------------------------------------------------
  {
    cout << "2D random hash test, doubles, coarse tolerance." << endl;
    const double xmin[2] = {-10.0, -10.0}, xmax[2] = {10.0, 10.0};
    const double tol = 1.5e-1;
    std::vector<double> coords(2*n);
    for (unsigned i = 0; i != n; ++i) {
      coords[2*i  ] = xmin[0] + random01()*(xmax[0] - xmin[0]);
      coords[2*i+1] = xmin[1] + random01()*(xmax[1] - xmin[1]);
    }
    hashSlingingSlasher<2, double>(xmin, xmax, tol, coords);
  }

  //----------------------------------------------------------------------
  // Hash a random set of points in 2D using a fine tolerance (ints).
  //----------------------------------------------------------------------
  {
    cout << "2D random hash test, int64_t, fine tolerance." << endl;
    const int64_t ixmax = (1 << 16);
    const int64_t xmin[2] = {-ixmax, -ixmax}, xmax[2] = {ixmax, ixmax};
    const int64_t tol = 1;
    std::vector<int64_t> coords(2*n);
    for (unsigned i = 0; i != n; ++i) {
      coords[2*i  ] = xmin[0] + int64_t(random01()*(xmax[0] - xmin[0]));
      coords[2*i+1] = xmin[1] + int64_t(random01()*(xmax[1] - xmin[1]));
    }
    hashSlingingSlasher<2, int64_t>(xmin, xmax, tol, coords);
  }

  //----------------------------------------------------------------------
  // Hash a random set of points in 2D using a coarse tolerance (ints).
  //----------------------------------------------------------------------
  {
    cout << "2D random hash test, int64_t, coarse tolerance." << endl;
    const int64_t ixmax = (1 << 16);
    const int64_t xmin[2] = {-ixmax, -ixmax}, xmax[2] = {ixmax, ixmax};
    const int64_t tol = 32;
    std::vector<int64_t> coords(2*n);
    for (unsigned i = 0; i != n; ++i) {
      coords[2*i  ] = xmin[0] + int64_t(random01()*(xmax[0] - xmin[0]));
      coords[2*i+1] = xmin[1] + int64_t(random01()*(xmax[1] - xmin[1]));
    }
    hashSlingingSlasher<2, int64_t>(xmin, xmax, tol, coords);
  }

  //----------------------------------------------------------------------
  // Hash a random set of points in 3D using a fine tolerance.
  //----------------------------------------------------------------------
  {
    cout << "3D random hash test, doubles, fine tolerance." << endl;
    const double xmin[3] = {-10.0, -10.0, -10.0}, xmax[3] = {10.0, 10.0, 10.0};
    const double tol = 1.5e-5;
    std::vector<double> coords(3*n);
    for (unsigned i = 0; i != n; ++i) {
      coords[3*i  ] = xmin[0] + random01()*(xmax[0] - xmin[0]);
      coords[3*i+1] = xmin[1] + random01()*(xmax[1] - xmin[1]);
      coords[3*i+2] = xmin[2] + random01()*(xmax[2] - xmin[2]);
    }
    hashSlingingSlasher<3, double>(xmin, xmax, tol, coords);
  }

  //----------------------------------------------------------------------
  // Hash a random set of points in 3D using a coarse tolerance.
  //----------------------------------------------------------------------
  {
    cout << "3D random hash test, doubles, coarse tolerance." << endl;
    const double xmin[3] = {-10.0, -10.0, -10.0}, xmax[3] = {10.0, 10.0, 10.0};
    const double tol = 1.5e-1;
    std::vector<double> coords(3*n);
    for (unsigned i = 0; i != n; ++i) {
      coords[3*i  ] = xmin[0] + random01()*(xmax[0] - xmin[0]);
      coords[3*i+1] = xmin[1] + random01()*(xmax[1] - xmin[1]);
      coords[3*i+2] = xmin[2] + random01()*(xmax[2] - xmin[2]);
    }
    hashSlingingSlasher<3, double>(xmin, xmax, tol, coords);
  }

  //----------------------------------------------------------------------
  // Hash a random set of points in 3D using a fine tolerance (ints).
  //----------------------------------------------------------------------
  {
    cout << "3D random hash test, int64_t, fine tolerance." << endl;
    const int64_t ixmax = (1 << 16);
    const int64_t xmin[3] = {-ixmax, -ixmax, -ixmax}, xmax[3] = {ixmax, ixmax, ixmax};
    const int64_t tol = 1;
    std::vector<int64_t> coords(3*n);
    for (unsigned i = 0; i != n; ++i) {
      coords[3*i  ] = xmin[0] + int64_t(random01()*(xmax[0] - xmin[0]));
      coords[3*i+1] = xmin[1] + int64_t(random01()*(xmax[1] - xmin[1]));
      coords[3*i+2] = xmin[2] + int64_t(random01()*(xmax[2] - xmin[2]));
    }
    hashSlingingSlasher<3, int64_t>(xmin, xmax, tol, coords);
  }

  //----------------------------------------------------------------------
  // Hash a random set of points in 3D using a coarse tolerance (ints).
  //----------------------------------------------------------------------
  {
    cout << "3D random hash test, int64_t, coarse tolerance." << endl;
    const int64_t ixmax = (1 << 16);
    const int64_t xmin[3] = {-ixmax, -ixmax, -ixmax}, xmax[3] = {ixmax, ixmax, ixmax};
    const int64_t tol = 33;
    std::vector<int64_t> coords(3*n);
    for (unsigned i = 0; i != n; ++i) {
      coords[3*i  ] = xmin[0] + int64_t(random01()*(xmax[0] - xmin[0]));
      coords[3*i+1] = xmin[1] + int64_t(random01()*(xmax[1] - xmin[1]));
      coords[3*i+2] = xmin[2] + int64_t(random01()*(xmax[2] - xmin[2]));
    }
    hashSlingingSlasher<3, int64_t>(xmin, xmax, tol, coords);
  }

  cout << "PASS" << endl;

#ifdef POLYTOPE_ENABLE_MPI
   MPI_Finalize();
#endif
  return 0;
}
