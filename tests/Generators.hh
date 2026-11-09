// This file must be brought in through polytope_test_utilities.hh, do not include it separately

#ifndef POLYTOPE_GENERATORS_HH
#define POLYTOPE_GENERATORS_HH

#include <algorithm>
#include <map>
#include <vector>
#include <random>

#include "polytope_boost_utilities.hh"
#include "Boundary2D.hh"
#include "Cube.hh"
#include "Communicator.hh"
#include "GeomUtils.hh"

using namespace std;

namespace polytope {
//------------------------------------------------------------------------
template<int Dimension>
class Generators {
public:
  // -------------- Public member variables and routines ---------------- //

  // Number of generators
  unsigned nPoints;
  vector<double> mPoints;
  Boundary2D* mBoundary = nullptr;
  Cube<double> mCube;

  //------------------------------------------------------------------------
  // Constructor, destructor
  //------------------------------------------------------------------------
  Generators(const Point<3, double>& min,
             const Point<3, double>& max) :
    mCube(min, max) {
  }

  Generators(Boundary2D& boundary) :
    nPoints(0),
    mPoints(0),
    mBoundary(&boundary) { };

  ~Generators() {};

  bool isInside(Point<Dimension, double>& pos) {
    if constexpr (Dimension == 2) {
      return boost::geometry::within(makeBGPoint(pos),
                                     mBoundary->mBGboundary);
    } else {
      return mCube.within(pos);
    }
  }

  //------------------------------------------------------------------------
  // Place random generators into spatial domain. Should be a uniform distribution
  //------------------------------------------------------------------------
  void randomPoints(unsigned nGenerators, unsigned seed = 0) {
    mPoints.clear();
    nPoints = nGenerators;
    srand(seed);
    // std::mt19937 gen(seed);
    // std::uniform_real_distribution<double> distrib(0., 1.);
    auto bHigh = mBoundary->maxBound;
    auto bLow = mBoundary->minBound;

    for (unsigned iter = 0; iter < nGenerators; ++iter ){
      Point<Dimension, double> pos;
      bool inside = false;
      while( !inside ) {
        for (unsigned n = 0; n < Dimension; ++n) {
          pos[n] = (bHigh[n]-bLow[n]) * random01() + bLow[n];
        }
        inside = isInside(pos);
      }
      for (auto n = 0; n < Dimension; ++n) {
        mPoints.push_back(pos[n]);
      }
    }
    POLY_CHECK( mPoints.size()/Dimension == nGenerators );
  }

  //------------------------------------------------------------------------
  // Place random generators into a 3D spatial domain.
  //------------------------------------------------------------------------
  void randomPoints(const Quantizer<Dimension>& Q,
                    unsigned nGenerators,
                    unsigned seed = 0) {
    mPoints.clear();
    nPoints = nGenerators;
    srand(seed);
    const auto bHigh = Q.domHiOrig();
    const auto bLow = Q.domLoOrig();

    for (unsigned iter = 0; iter < nGenerators; ++iter) {
      Point<Dimension, double> pos;
      bool inside = false;
      while (!inside) {
        for (unsigned n = 0; n < Dimension; ++n) {
          pos[n] = (bHigh[n] - bLow[n])*random01() + bLow[n];
        }
        inside = isInside(pos);
      }
      for (unsigned n = 0; n < Dimension; ++n) mPoints.push_back(pos[n]);
    }
    POLY_CHECK(mPoints.size()/Dimension == nGenerators);
  }

  //------------------------------------------------------------------------
  // Place random generators into spatial domain using a normal distribution
  // Focuses more points near the center of the domain
  //------------------------------------------------------------------------
  void randomNormalPoints(unsigned nGenerators, unsigned seed = 0) {
    mPoints.clear();
    nPoints = nGenerators;
    double mean = 0.5;
    double sigma = 1./6.17737; // Provides a 99.9% confidence interval to fall between 0 and 1
    std::mt19937 gen(seed);
    std::normal_distribution<double> distrib(mean, sigma);
    auto bHigh = mBoundary->maxBound;
    auto bLow = mBoundary->minBound;

    for (unsigned iter = 0; iter < nGenerators; ++iter ){
      Point<Dimension, double> pos;
      bool inside = false;
      while( !inside ) {
        for (unsigned n = 0; n < Dimension; ++n) {
          pos[n] = (bHigh[n]-bLow[n]) * distrib(gen) + bLow[n];
        }
        inside = isInside(pos);
      }
      for (auto n = 0; n < Dimension; ++n) {
        mPoints.push_back(pos[n]);
      }
    }
    POLY_CHECK( mPoints.size()/Dimension == nGenerators );
  }

  //------------------------------------------------------------------------
  // Place Cartesian points of constant mesh spacing
  //------------------------------------------------------------------------
  void cartesianPoints(vector<unsigned> nCellsPerDimension) {
    mPoints.clear();
    POLY_CHECK( nCellsPerDimension.size() == Dimension );

    if constexpr (Dimension == 2) {
      cartesian2D(nCellsPerDimension[0], nCellsPerDimension[1]);
    } else {
      POLY_CHECK2(false, "Cartesian generation is only supported for 2D boundaries");
    }
    nPoints = mPoints.size()/Dimension;
  }

  //------------------------------------------------------------------------
  // add a point to the generator set
  //------------------------------------------------------------------------
  void addGenerator(double* pos) {
    std::vector<double> point;
    for (unsigned n=0; n<Dimension; ++n ) point.push_back( pos[n] );
    bool inside = boost::geometry::within( makeBGPoint(point), mBoundary->mBGboundary );
    POLY_CHECK( inside );
    mPoints.insert( mPoints.end(), point.begin(), point.end() );
  }

  //------------------------------------------------------------------------
  // for 2D problems
  //------------------------------------------------------------------------
  void cartesian2D(const unsigned nx, const unsigned ny) {
    auto bHigh = mBoundary->maxBound;
    auto bLow = mBoundary->minBound;
    double x, y;
    double dx = (bHigh[0] - bLow[0]) / nx;
    double dy = (bHigh[1] - bLow[1]) / ny;
    Point2<double> pos;
    for (unsigned iy = 0; iy != ny; ++iy) {
      y = bLow[1] + (iy + 0.5)*dy;
      for (unsigned ix = 0; ix != nx; ++ix) {
        x = bLow[0] + (ix + 0.5)*dx;
        pos[0] = x;
        pos[1] = y;
        if( boost::geometry::within( makeBGPoint(pos), mBoundary->mBGboundary) ){
          mPoints.push_back( x );
          mPoints.push_back( y );
        }
      }
    }
    nPoints = mPoints.size()/Dimension;
  }

  //------------------------------------------------------------------------
  // perturb the locations of the generators by +/- epsilon/2 in each direction
  //------------------------------------------------------------------------
  void perturb(double epsilon) {
    for (unsigned i = 0; i < mPoints.size()/Dimension; ++i){
      for (unsigned n = 0; n < Dimension; ++n){
        mPoints[Dimension*i+n] += epsilon*(random01() - 0.5 );
      }
    }
  }

  //------------------------------------------------------------------------
  // Create a Boost.Geometry point from a std::vector of data depending on
  // the dimension of the problem. There really should be an easier way
  //------------------------------------------------------------------------
  BGPoint<double, Dimension> makeBGPoint(Point<Dimension, double> pointIn) {
    if constexpr (Dimension == 2) {
      return BGPoint<double, 2>(pointIn[0], pointIn[1]);
    } else if constexpr (Dimension == 3) {
      return BGPoint<double, 3>(pointIn[0], pointIn[1], pointIn[2]);
    }
  }
};
}
#endif
