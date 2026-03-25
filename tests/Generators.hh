// This file must be brought in through polytope_test_utilities.hh, do not include it separately

#ifndef POLYTOPE_GENERATORS_HH
#define POLYTOPE_GENERATORS_HH

#include "Boundary2D.hh"

using namespace std;

namespace polytope {
//------------------------------------------------------------------------
template<int Dimension, typename RealType>
class Generators {
public:
   // -------------- Public member variables and routines ---------------- //

   // Number of generators
   unsigned nPoints;
   vector<RealType> mPoints;
   Boundary2D<RealType>& mBoundary;   
   
   //------------------------------------------------------------------------
   // Constructor, destructor
   //------------------------------------------------------------------------
   Generators(Boundary2D<RealType>& boundary):
      nPoints(0),
      mPoints(0),
      mBoundary(boundary) {};

   ~Generators() {};


   //------------------------------------------------------------------------
   // Place random generators into spatial domain
   //------------------------------------------------------------------------
   void randomPoints(unsigned nGenerators)
   {
      mPoints.clear();
      nPoints = nGenerators;

      mBoundary.getBoundingBox();

      for (unsigned iter = 0; iter < nGenerators; ++iter ){
         std::vector<RealType> pos(Dimension,0);
         bool inside = false;
         while( !inside ){
            for (unsigned n = 0; n < Dimension; ++n){
               pos[n] = (mBoundary.mHigh[n]-mBoundary.mLow[n]) * 
                  RealType(::random())/RAND_MAX + mBoundary.mLow[n];
            }
            inside = boost::geometry::within( makeBGPoint(pos),
                                              mBoundary.mBGboundary );
         }
         mPoints.insert( mPoints.end(), pos.begin(), pos.end() );
      }
      POLY_ASSERT( mPoints.size()/Dimension == nGenerators );
   }

   
   //------------------------------------------------------------------------
   // Place Cartesian points of constant mesh spacing
   //------------------------------------------------------------------------
   void cartesianPoints(vector<unsigned> nCellsPerDimension)
   {
      mPoints.clear();
      POLY_ASSERT( nCellsPerDimension.size() == Dimension );

      mBoundary.getBoundingBox();

      if( Dimension == 2 ){
         cartesian2D( nCellsPerDimension[0], nCellsPerDimension[1] );
      }else{
         cartesian3D( nCellsPerDimension[0], nCellsPerDimension[1], nCellsPerDimension[2] );
      }
      nPoints = mPoints.size()/Dimension;
   }
    

   //------------------------------------------------------------------------
   // Place radial generators about the center specified in Boundary2D
   //------------------------------------------------------------------------
   void radialPoints(const unsigned nr)
   {
      mPoints.clear();
      POLY_ASSERT( Dimension == 2 );
      RealType maxDistance;
      mBoundary.getBoundingRadius( maxDistance );
      POLY_ASSERT( maxDistance > 0 );
      
      RealType dRadius = maxDistance/nr;
      for( unsigned i = 0; i != nr; ++i ){
         RealType rad = (i+0.5)*dRadius;

         // This is supposed to befloor(2*pi*i), however 2*floor(pi)*i=6*i 
         // is found to work better
         unsigned nArcs = 6*i;
         std::vector<RealType> pos(2,0);
         for( unsigned j = 0; j != nArcs; ++j ){
            RealType theta = 2*M_PI*j/nArcs;
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
   void addGenerator(RealType* pos)
   {
      std::vector<RealType> point;
      for (unsigned n=0; n<Dimension; ++n ) point.push_back( pos[n] );
      bool inside = boost::geometry::within( makeBGPoint(point), mBoundary.mBGboundary );
      POLY_ASSERT( inside );
      mPoints.insert( mPoints.end(), point.begin(), point.end() );
   }

   //------------------------------------------------------------------------
   // for 2D problems
   //------------------------------------------------------------------------
   void cartesian2D(const unsigned nx, const unsigned ny)
   {
      RealType x, y;
      RealType dx = (mBoundary.mHigh[0] - mBoundary.mLow[0]) / nx;
      RealType dy = (mBoundary.mHigh[1] - mBoundary.mLow[1]) / ny;
      std::vector<RealType> pos(Dimension,0);
      for (unsigned iy = 0; iy != ny; ++iy){
         y = mBoundary.mLow[1] + (iy + 0.5)*dy;
         for (unsigned ix = 0; ix != nx; ++ix){
            x = mBoundary.mLow[0] + (ix + 0.5)*dx;
            pos[0] = x;   pos[1] = y;
            if( boost::geometry::within( makeBGPoint(pos), mBoundary.mBGboundary) ){
               mPoints.push_back( x );
               mPoints.push_back( y );
            }
         }
      }
   }

   //------------------------------------------------------------------------
   // for 3D problems
   //------------------------------------------------------------------------
   void cartesian3D(const unsigned nx, const unsigned ny, const unsigned nz)
   {
      RealType x, y, z;
      RealType dx = (mBoundary.mHigh[0] - mBoundary.mLow[0]) / nx;
      RealType dy = (mBoundary.mHigh[1] - mBoundary.mLow[1]) / ny;
      RealType dz = (mBoundary.mHigh[2] - mBoundary.mLow[2]) / nz;
      for (unsigned iz = 0; iz != nz; ++iz){
         z = mBoundary.mLow[2] + (iz + 0.5)*dz;
         for (unsigned iy = 0; iy != ny; ++iy){
            y = mBoundary.mLow[1] + (iy + 0.5)*dy;
            for (unsigned ix = 0; ix != nx; ++ix){
               x = mBoundary.mLow[0] + (ix + 0.5)*dx;
               mPoints.push_back( x );
               mPoints.push_back( y );
               mPoints.push_back( z );
            }
         }
      }
   }

   //------------------------------------------------------------------------
   // perturb the locations of the generators by +/- epsilon/2 in each direction
   //------------------------------------------------------------------------
   void perturb(RealType epsilon)
   {
      for (unsigned i = 0; i < mPoints.size()/Dimension; ++i){
         for (unsigned n = 0; n < Dimension; ++n){
            mPoints[Dimension*i+n] += epsilon*( RealType(::random())/RAND_MAX - 0.5 );
         }
      }
   }

  //------------------------------------------------------------------------
  // Create a Boost.Geometry point from a std::vector of data depending on
  // the dimension of the problem. There really should be an easier way
  //------------------------------------------------------------------------
  BGPoint<RealType, Dimension> makeBGPoint(std::vector<RealType> pointIn) {
    if constexpr (Dimension == 2) {
        return BGPoint<RealType, 2>(pointIn[0], pointIn[1]);
      } else if constexpr (Dimension == 3) {
        return BGPoint<RealType, 3>(pointIn[0], pointIn[1], pointIn[2]);
      }
  }
};
}
#endif
