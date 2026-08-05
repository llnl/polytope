// This file must be brought in through polytope_test_utilities.hh, do not include it separately

#ifndef POLYTOPE_GENERATORS_HH
#define POLYTOPE_GENERATORS_HH

#include "polytope_boost_utilities.hh"
#include "Boundary2D.hh"
#include "Cube.hh"
#include "Communicator.hh"
#include "GeomUtils.hh"

#include <algorithm>
#include <map>
#include <vector>

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
  Boundary2D mBoundary;
  Cube<double> mCube;

  //------------------------------------------------------------------------
  // Constructor, destructor
  //------------------------------------------------------------------------
  Generators(const Point<3, double>& min,
             const Point<3, double>& max) :
    mCube(min, max) {
    auto& Q = Quantizer<3>::instance();
    Q.init(min, max);
  }

  Generators(Boundary2D& boundary) :
    nPoints(0),
    mPoints(0),
    mBoundary(boundary) {};

  ~Generators() {};

  bool isInside(Point<Dimension, double>& pos) {
    if constexpr (Dimension == 2) {
      return boost::geometry::within(makeBGPoint(pos),
                                     mBoundary.mBGboundary);
    } else {
      return mCube.within(pos);
    }
  }

  //------------------------------------------------------------------------
  // Place random generators into spatial domain
  //------------------------------------------------------------------------
  void randomPoints(unsigned nGenerators, unsigned seed = 0) {
    mPoints.clear();
    nPoints = nGenerators;
    srand(seed);
    // std::mt19937 gen(seed);
    // std::uniform_real_distribution<double> distrib(0., 1.);
    auto& Q = Quantizer<Dimension>::instance();
    auto bHigh = Q.m_xhi;
    auto bLow = Q.m_xlo;

    for (unsigned iter = 0; iter < nGenerators; ++iter ){
      Point<Dimension, double> pos;
      bool inside = false;
      while( !inside ) {
        for (unsigned n = 0; n < Dimension; ++n){
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
  // Place Cartesian points of constant mesh spacing
  //------------------------------------------------------------------------
  void cartesianPoints(vector<unsigned> nCellsPerDimension) {
    mPoints.clear();
    POLY_CHECK( nCellsPerDimension.size() == Dimension );

    if (Dimension == 2){
      cartesian2D( nCellsPerDimension[0], nCellsPerDimension[1] );
    } else {
      cartesian3D( nCellsPerDimension[0], nCellsPerDimension[1], nCellsPerDimension[2] );
    }
    nPoints = mPoints.size()/Dimension;
  }

  //------------------------------------------------------------------------
  // Place radial generators about the center specified in Boundary2D
  //------------------------------------------------------------------------
  void radialPoints(const unsigned nr) {
    mPoints.clear();
    POLY_CHECK( Dimension == 2 );
    double maxDistance;
    mBoundary.getBoundingRadius( maxDistance );
    POLY_CHECK( maxDistance > 0 );

    double dRadius = maxDistance/nr;
    for( unsigned i = 0; i != nr; ++i ) {
      double rad = (i+0.5)*dRadius;

      // This is supposed to befloor(2*pi*i), however 2*floor(pi)*i=6*i
      // is found to work better
      unsigned nArcs = 6*i;
      std::vector<double> pos(2,0);
      for( unsigned j = 0; j != nArcs; ++j ){
        double theta = 2*M_PI*j/nArcs;
        pos[0] = mBoundary.mCenter[0] + rad*cos(theta);
        pos[1] = mBoundary.mCenter[1] + rad*sin(theta);
        if( boost::geometry::within( makeBGPoint(pos), mBoundary.mBGboundary) ){
          mPoints.push_back( pos[0] );
          mPoints.push_back( pos[1] );
        }
      }
    }
    nPoints = mPoints.size()/Dimension;
  }

  //------------------------------------------------------------------------
  // add a point to the generator set
  //------------------------------------------------------------------------
  void addGenerator(double* pos) {
    std::vector<double> point;
    for (unsigned n=0; n<Dimension; ++n ) point.push_back( pos[n] );
    bool inside = boost::geometry::within( makeBGPoint(point), mBoundary.mBGboundary );
    POLY_CHECK( inside );
    mPoints.insert( mPoints.end(), point.begin(), point.end() );
  }

  //------------------------------------------------------------------------
  // for 2D problems
  //------------------------------------------------------------------------
  void cartesian2D(const unsigned nx, const unsigned ny) {
    auto& Q = Quantizer<2>::instance();
    auto bHigh = Q.m_xhi;
    auto bLow = Q.m_xlo;
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
        if( boost::geometry::within( makeBGPoint(pos), mBoundary.mBGboundary) ){
          mPoints.push_back( x );
          mPoints.push_back( y );
        }
      }
    }
    nPoints = mPoints.size()/Dimension;
  }

  //------------------------------------------------------------------------
  // for 3D problems
  //------------------------------------------------------------------------
  void cartesian3D(const unsigned nx, const unsigned ny, const unsigned nz) {
    auto& Q = Quantizer<3>::instance();
    auto bHigh = Q.m_xhi;
    auto bLow = Q.m_xlo;
    double x, y, z;
    double dx = (bHigh[0] - bLow[0]) / nx;
    double dy = (bHigh[1] - bLow[1]) / ny;
    double dz = (bHigh[2] - bLow[2]) / nz;
    for (unsigned iz = 0; iz != nz; ++iz) {
      z = bLow[2] + (iz + 0.5)*dz;
      for (unsigned iy = 0; iy != ny; ++iy) {
        y = bLow[1] + (iy + 0.5)*dy;
        for (unsigned ix = 0; ix != nx; ++ix) {
          x = bLow[0] + (ix + 0.5)*dx;
          mPoints.push_back( x );
          mPoints.push_back( y );
          mPoints.push_back( z );
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

  //------------------------------------------------------------------------
  // Parallel utilities
  //------------------------------------------------------------------------

  //------------------------------------------------------------------------
  // Assigns a random index from 0 to nPoints for each proc between 0 and maxProc.
  // maxProc can cap which processors are provided points. Defaults to all procs
  //------------------------------------------------------------------------
  std::vector<unsigned> assignRandomPointToRank(unsigned seed = 0, int maxProc = -1) {
    srand(seed);
    nPoints = mPoints.size()/Dimension;
    // Figure out parallel configuration
    int numProcs = Communicator::getNProcs();

    if (maxProc >= 0) {
      numProcs = maxProc;
    }
    numProcs = std::min(int(nPoints), numProcs);
    std::vector<unsigned> pointIndices(nPoints);
    for (unsigned i = 0; i < nPoints; ++i) pointIndices[i] = i;

    // For each proc, grab a random index from 0 to nPoints.
    // That proc gets that point as it's centralized location.
    std::vector<unsigned> procIndex(numProcs);
    for (unsigned iproc = 0; iproc < static_cast<unsigned>(numProcs); ++iproc) {
      const auto remaining = nPoints - iproc;
      auto offset = static_cast<unsigned>(random01()*remaining);
      offset = std::min(offset, remaining - 1);
      std::swap(pointIndices[iproc], pointIndices[iproc + offset]);
      procIndex[iproc] = pointIndices[iproc];
    }
    return procIndex;
  }

  //------------------------------------------------------------------------
  // Assign points based on processor by removing non-local points
  // Returns a map from each original point to its assigned rank
  //------------------------------------------------------------------------
  std::map<Point<Dimension, double>, int> distributePointsAmongRanks(const std::vector<unsigned>& procIndex) {
    auto& Q = Quantizer<Dimension>::instance();
    std::map<Point<Dimension, double>, int> finalRanks;
    // Figure out parallel configuration
    int rank = Communicator::getRank();
    int numProcs = procIndex.size();

    std::vector<Point<Dimension, double>> procPoint(numProcs);
    for (int iproc = 0; iproc < numProcs; ++iproc) {
      const auto pin = procIndex[iproc];
      for (int d = 0; d < Dimension; ++d) {
        procPoint[iproc][d] = mPoints[Dimension*pin+d];
      }
    }

    std::vector<double> newPoints;
    newPoints.reserve(mPoints.size());
    for (int i = 0; i < nPoints; ++i) {
      int owner = 0;
      Point<Dimension, double> point;
      for (int d = 0; d < Dimension; ++d) {
        point[d] = mPoints[Dimension*i+d];
      }
      auto diff = point - procPoint[0];
      double minDist = magnitude<Dimension>(diff);
      for (int ip = 1; ip < numProcs; ++ip) {
        diff = point - procPoint[ip];
        double dis = magnitude<Dimension>(diff);
        if (dis < minDist) {
          owner = ip;
          minDist = dis;
        }
      }
      // Do a round-trip quant/dequant on the point to ensure it matches
      // the point in the tessellated mesh.
      Point<Dimension, double> qpoint = Q.dequantize(Q.quantize(point));
      finalRanks[qpoint] = owner;
      if (owner == rank) {
        for (int d = 0; d < Dimension; ++d) {
          newPoints.push_back(point[d]);
        }
      }
    }
    mPoints = std::move(newPoints);
    nPoints = mPoints.size()/Dimension;
    return finalRanks;
  }
};
}
#endif
