#ifndef POLYTOPE_BOUNDARY2D_HH
#define POLYTOPE_BOUNDARY2D_HH

//#include "within.hh"
#include <iostream>
#include "polytope_boost_utilities.hh"
#include "QuantPLC.hh"
#include "Shapes.hh"

//------------------------------------------------------------------------------
// Define our own local random number generator wrapping the standard srand &
// rand methods.
//------------------------------------------------------------------------------
inline
double random01() {
  return double(rand())/RAND_MAX;
}

using namespace std;
namespace polytope {
//------------------------------------------------------------------------
template<typename RealType>
class Boundary2D {
public:
  // -------------- Public member variables and routines ---------------- //
  
  // Piecewise linear construct to define the boundary facets + holes
  PLC<2> mPLC;
  // Vector of generators to define the boundary
  std::vector<RealType> mPLCpoints;
  double mDiff = 0.5;
  // Ranges of bounding box
  RealType mCenter[3] = {0.0, 0.0, 0.0};
  RealType mArea;

  BGPolygon<RealType, 2> mBGboundary;
  
  // Define enum to keep track fo the type of boundary called for
  enum BoundaryType{
    square             = 0,
    circle             = 1,
    donut              = 2,
    mwithholes         = 3,
    funkystar          = 4,
    circlewithstarhole = 5,
    cardioid           = 6,
    trogdor            = 7,
    starwithhole       = 8,
    trogdor2           = 9
  };
  
  // Boundary type
  mutable BoundaryType mType;

  double m_pad = 0.1;
  Quantizer<2> mQ;
  QuantPLC<2> mQPLC;
  
  //------------------------------------------------------------------------
  // Constructor, destructor
  //------------------------------------------------------------------------
  Boundary2D():
    mType(square){
    this->clear();
    this->setDefaultBoundary(mType);
  }
  
  ~Boundary2D() {};
    
  void clear() {
    mPLC.clear();
    mPLCpoints.clear();
    boost::geometry::clear(mBGboundary);
  }
  
  void finalize() {
    mQ.init(mPLCpoints, -1., m_pad);
    mQPLC.init(mPLC, mQ, mPLCpoints);
    boostMyBoundary();
    mArea = boost::geometry::area(mBGboundary);
  }
  
  //------------------------------------------------------------------------
  // setCustomBoundary
  //------------------------------------------------------------------------
  void setCustomBoundary(const unsigned numVertices,
                         const RealType* vertices) {
    this->clear();
    mPLCpoints.resize(2*numVertices);
    mPLC.facets.resize(numVertices, vector<int>(2));
    for (unsigned i = 0; i != numVertices; ++i) {
      mPLCpoints[2*i  ] = vertices[2*i  ];
      mPLCpoints[2*i+1] = vertices[2*i+1];
      mPLC.facets[i][0] = i;
      mPLC.facets[i][1] = (i+1) % numVertices;
    }
    
    this->finalize();
  }

  //------------------------------------------------------------------------
  // setDefaultBoundary
  //------------------------------------------------------------------------
  void setDefaultBoundary(const int bType)
  {

    switch(bType){
    case square:
      this->setUnitSquare();
      break;
    case circle:
      this->setUnitCircle();
      break;
    case donut:
      this->setDonut();
      break;
    case mwithholes:
      this->setMWithHoles();
      break;
    case funkystar:
      this->setFunkyStar();
      break;
    case circlewithstarhole:
      this->setCircleWithStarHole();
      break;
    case cardioid:
      this->setCardioid();
      break;
    case trogdor:
      this->setTrogdor();
      break;
    case starwithhole:
      this->setStarWithHole();
      break;
    case trogdor2:
      this->setTrogdor2();
      break;
    }
   }

  //------------------------------------------------------------------------
  // setUnitSquare
  //------------------------------------------------------------------------
  void setUnitSquare() {
    this->clear();

    const RealType x1 = mCenter[0] - mDiff;
    const RealType x2 = mCenter[0] + mDiff;
    const RealType y1 = mCenter[1] - mDiff;
    const RealType y2 = mCenter[1] + mDiff;
    Point2<RealType> low(x1, y1);
    Point2<RealType> hi(x2, y2);
    std::vector<Point2<RealType>> points = shapes::createSquarePoints(low, hi);
    mPLCpoints = flattenCoords(points);
    mPLC.facets = shapes::createSquareFaces();
    mType = square;
    this->finalize();
  }

  //------------------------------------------------------------------------
  // setUnitCircle
  //------------------------------------------------------------------------
  void setUnitCircle() {
    this->clear();
    // Boundary generators.
    unsigned Nb = 90; // 4-degree resolution.
    for (unsigned b = 0; b < Nb; ++b) {
      RealType theta = 2.0*M_PI*b/(Nb+1);
      RealType x = mCenter[0] + cos(theta);
      RealType y = mCenter[1] + sin(theta);
      mPLCpoints.push_back(x);
      mPLCpoints.push_back(y);
    }
    
    // Facets.
    mPLC.facets.resize(Nb); 
    for (unsigned f = 0; f < Nb; ++f) {
      unsigned fBegin =  mPLCpoints.size()/2 - Nb + f;
      unsigned fEnd   = (mPLCpoints.size()/2 - Nb + f + 1) % Nb;
      mPLC.facets[f].resize(2);
      mPLC.facets[f][0] = fBegin;
      mPLC.facets[f][1] = fEnd;
    }
    mType = circle;
    this->finalize();
  }
  
  //------------------------------------------------------------------------
  // setDonut
  // Unit circle with a hole of prescribed radius
  //------------------------------------------------------------------------
  void setDonut( RealType innerRadius=0.25 ) {
    POLY_CHECK2( innerRadius > 0, "Must provide a positive inner radius" );
    POLY_CHECK2( innerRadius < 1, "Inner radius may not exceed outer (unit) radius" );
    
    // The outer circle
    this->clear();
    this->setUnitCircle();
    
    // Inner circle.
    unsigned Nb = mPLCpoints.size()/2;
    unsigned Nh = Nb;
    for (unsigned b = 0; b < Nh; ++b) {
      RealType theta = 2.0*M_PI*(1.0 - RealType(b)/RealType(Nb+1));
      RealType x = mCenter[0] + innerRadius*cos(theta);
      RealType y = mCenter[1] + innerRadius*sin(theta);
      mPLCpoints.push_back(x);
      mPLCpoints.push_back(y);
    }
    
    // Facets on the inner circle
    mPLC.holes = vector< vector< vector<int> > >(1);
    mPLC.holes[0].resize(Nh);
    for (unsigned f = 0; f < Nh; ++f) {
      unsigned fBegin = mPLCpoints.size()/2 - Nh + f;
      unsigned fEnd   = mPLCpoints.size()/2 - Nh + (f + 1) % Nh;
      mPLC.holes[0][f].resize(2);
      mPLC.holes[0][f][0] = fBegin;
      mPLC.holes[0][f][1] = fEnd;
    }
    mType = donut;
    this->finalize();
  }
  
  //------------------------------------------------------------------------
  // setMWithHoles
  // M-shaped domain with two square holes
  //------------------------------------------------------------------------
  void setMWithHoles() {
    this->clear();
    // Outer boundary of the M-shape
    mPLCpoints.push_back(0.0); mPLCpoints.push_back(0.0);
    mPLCpoints.push_back(2.0); mPLCpoints.push_back(0.0);
    mPLCpoints.push_back(2.0); mPLCpoints.push_back(2.0);
    mPLCpoints.push_back(1.0); mPLCpoints.push_back(1.0);
    mPLCpoints.push_back(0.0); mPLCpoints.push_back(2.0);
    
    int nSides = mPLCpoints.size()/2;
    mPLC.facets.resize( nSides, vector<int>(2) );
    for (unsigned i = 0; i != nSides; ++i ) {
      mPLC.facets[i][0] = i;
      mPLC.facets[i][1] = (i+1) % nSides;
    }
    
    int nHoles = 2;
    
    // Square hole #1
    mPLCpoints.push_back(0.25); mPLCpoints.push_back(0.25);
    mPLCpoints.push_back(0.25); mPLCpoints.push_back(0.75);
    mPLCpoints.push_back(0.75); mPLCpoints.push_back(0.75);
    mPLCpoints.push_back(0.75); mPLCpoints.push_back(0.25);
    
    // Square hole #2
    mPLCpoints.push_back(1.25); mPLCpoints.push_back(0.25);
    mPLCpoints.push_back(1.25); mPLCpoints.push_back(0.75);
    mPLCpoints.push_back(1.75); mPLCpoints.push_back(0.75);
    mPLCpoints.push_back(1.75); mPLCpoints.push_back(0.25);
    
    mPLC.holes.resize(nHoles, vector<vector<int> >(4, vector<int>(nHoles)));
    for (unsigned i = 0; i != 4; ++i) {
      mPLC.holes[0][i][0] = nSides + i;
      mPLC.holes[0][i][1] = nSides + ((i+1) % 4);
      mPLC.holes[1][i][0] = nSides + 4 + i;
      mPLC.holes[1][i][1] = nSides + 4 + ((i+1) % 4);
    }
    mType = mwithholes;
    this->finalize();
  }
  
  //------------------------------------------------------------------------
  // setFunkyStar
  // Star-shaped(-ish) region from Misha Shashkov's Voronoi test suite
  //------------------------------------------------------------------------
  void setFunkyStar() {
    this->clear();
    // Get the boundary points, organized counterclockwise
    int Nsides = 10;
    for (int i = Nsides-1; i >= 0; --i ) {
      RealType rad = 2.0 + 0.5*sin(12.0* M_PI*RealType(i)/(RealType(Nsides)-1.0));
      RealType x = rad * sin(M_PI*RealType(i)/(RealType(Nsides)-1.0));
      RealType y = rad * cos(M_PI*RealType(i)/(RealType(Nsides)-1.0));
      mPLCpoints.push_back( x );
      mPLCpoints.push_back( y );
    }
    
    // Connect the boundary facets
    mPLC.facets.resize(Nsides);
    for (unsigned f = 0; f < Nsides; ++f ) {
      mPLC.facets[f].resize(2);
      mPLC.facets[f][0] = f;
      mPLC.facets[f][1] = (f+1) % Nsides;
    }
    mType = funkystar;
    this->finalize();
  }

  //------------------------------------------------------------------------
  // circleWithStarHole
  // Unit circle with a hole shaped like a regular n-pointed star
  //------------------------------------------------------------------------
  void setCircleWithStarHole( int nPoints = 5 ) {
    this->clear();
    // The outer boundary
    this->setUnitCircle();
    
    RealType theta0 = 2*M_PI/nPoints;
    RealType outerRadius = 0.75;
    RealType innerRadius = outerRadius*( sin(theta0/4.0) / sin(3*theta0/4.0) );
    
    RealType theta;
    for (unsigned p = 0; p < nPoints; ++p ) {
      // For the pointy bits of the star
      theta = M_PI/2 - p*theta0;
      mPLCpoints.push_back( mCenter[0] + outerRadius*cos(theta) );
      mPLCpoints.push_back( mCenter[1] + outerRadius*sin(theta) );
      
      // For the concave bits of the star
      theta = M_PI/2 - p*theta0 - theta0/2.0;
      mPLCpoints.push_back( mCenter[0] + innerRadius*cos(theta) );
      mPLCpoints.push_back( mCenter[1] + innerRadius*sin(theta) );
    }
    
    // Facets on the inner circle
    mPLC.holes = vector< vector< vector<int> > >(1);      
    mPLC.holes[0].resize(2*nPoints);
    for (unsigned f = 0; f < 2*nPoints; ++f) {
      unsigned fBegin = mPLCpoints.size()/2 - 2*nPoints + f;
      unsigned fEnd   = mPLCpoints.size()/2 - 2*nPoints + ((f + 1) % (2*nPoints));
      mPLC.holes[0][f].resize(2);
      mPLC.holes[0][f][0] = fBegin;
      mPLC.holes[0][f][1] = fEnd;
    }
    
    mType = circlewithstarhole;
    this->finalize();
  }
   
  //------------------------------------------------------------------------
  // cardioid
  // Cardioid with parametric equations
  //   x(t) = z*cos(t) + cos(2*t)
  //   y(t) = z*sin(t) + sin(2*t)
  // for t in [0,2*pi] and coefficient z
  //
  // NOTE:   z=2 : defines the unit cardioid with cusp pointing inward
  //       0<z<2 : cusp pushes inward and creates a hole
  //         z>2 : cusp less pronounced, cardioid more circular, non-unit
  //------------------------------------------------------------------------
  void setCardioid( RealType z = 2 ) {
    POLY_CHECK2( z > 0, "Must provide a positive coefficient for the cardioid" );
    this->clear();
    
    // Boundary generators.
    unsigned Nb = 90; // 4-degree resolution.
    for (unsigned b = 0; b < Nb; ++b) {
      RealType theta = 2.0*M_PI*b/(Nb+1);
      RealType x = z*cos(theta) - cos(2*theta);
      RealType y = z*sin(theta) - sin(2*theta);
      mPLCpoints.push_back(x);
      mPLCpoints.push_back(y);
    }
    
    // Facets.
    mPLC.facets.resize(Nb); 
    for (unsigned f = 0; f < Nb; ++f) {
      unsigned fBegin =  mPLCpoints.size()/2 - Nb + f;
      unsigned fEnd   = (mPLCpoints.size()/2 - Nb + f + 1) % Nb;
      mPLC.facets[f].resize(2);
      mPLC.facets[f][0] = fBegin;
      mPLC.facets[f][1] = fEnd;
    }
    mType = cardioid;
    this->finalize();
  }
   
  //------------------------------------------------------------------------
  // setTrogdor
  // No explanation necessary
  //------------------------------------------------------------------------
  void setTrogdor() {
    this->clear();
    const unsigned nSides = 30;
    const RealType points[60] = {2.0, 9.0, 4.0, 8.9, 5.0, 9.2, 6.5, 8.8, 
				 7.0, 8.0, 6.5, 7.0, 5.0, 6.3, 4.0, 5.5, 
				 3.7, 4.6, 4.0, 3.5, 5.5, 2.8, 6.8, 3.4,
				 5.5, 2.5, 4.0, 2.6, 3.0, 3.0, 2.5, 4.0,
				 2.4, 4.5, 2.5, 5.3, 3.0, 5.9, 4.5, 7.0,
				 4.9, 7.3, 5.1, 7.7, 5.0, 8.0, 4.5, 8.0,
				 3.5, 7.5, 2.0, 7.0, 2.0, 7.4, 3.3, 8.0,
				 1.75, 8.0, 1.75, 9.2};
    
    mPLC.facets.resize( nSides, vector<int>(2) );
    for (unsigned i = 0; i != nSides; ++i){
      unsigned j = nSides - i - 1;
      mPLCpoints.push_back(points[2*j  ]);
      mPLCpoints.push_back(points[2*j+1]);
      mPLC.facets[i][0] = i;
      mPLC.facets[i][1] = (i+1) % nSides;
    }
    
    mType = trogdor;
    this->finalize();
  }

  //------------------------------------------------------------------------
  // setStarWithHole
  // 5-pt star with hole in center from Misha Shashkov's Voronoi test suite
  //------------------------------------------------------------------------
  void setStarWithHole() {
    this->clear();
    const unsigned nPoints = 5;
    const RealType theta0 = 2*M_PI/nPoints;
    const RealType outerRadius = 1.0;
    const RealType innerRadius = outerRadius*( sin(theta0/4.0) / sin(3*theta0/4.0) );
    
    RealType theta;
    for (unsigned p = 0; p < nPoints; ++p ) {
      // For the pointy bits of the star
      theta = M_PI/2 + p*theta0;
      mPLCpoints.push_back( mCenter[0] + outerRadius*cos(theta) );
      mPLCpoints.push_back( mCenter[1] + outerRadius*sin(theta) );
      
      // For the concave bits of the star
      theta = M_PI/2 + p*theta0 + theta0/2.0;
      mPLCpoints.push_back( mCenter[0] + innerRadius*cos(theta) );
      mPLCpoints.push_back( mCenter[1] + innerRadius*sin(theta) );
    }
    
    // Facets on the inner circle
    mPLC.facets.resize( 2*nPoints, vector<int>(2) );
    for (unsigned f = 0; f < 2*nPoints; ++f) {
      unsigned fBegin = mPLCpoints.size()/2 - 2*nPoints + f;
      unsigned fEnd   = mPLCpoints.size()/2 - 2*nPoints + ((f + 1) % (2*nPoints));
      mPLC.facets[f][0] = fBegin;
      mPLC.facets[f][1] = fEnd;
    }
    
    // The points that define the inner hole
    const unsigned nHolePoints = 4;
    const RealType holePoints[8] = {0.05, -0.05,
				    0.10,  0.10,
				    0.20, -0.30,
				   -0.25, -0.15};
    for (unsigned p = 0; p < nHolePoints; ++p){
      mPLCpoints.push_back( holePoints[2*p  ] );
      mPLCpoints.push_back( holePoints[2*p+1] );
    }
    
    // Facets on the inner circle
    mPLC.holes = vector< vector< vector<int> > >(1);      
    mPLC.holes[0].resize(nHolePoints);
    for (unsigned f = 0; f < nHolePoints; ++f) {
      unsigned fBegin = mPLCpoints.size()/2 - nHolePoints + f;
      unsigned fEnd   = mPLCpoints.size()/2 - nHolePoints + ((f + 1) % nHolePoints);
      mPLC.holes[0][f].resize(2);
      mPLC.holes[0][f][0] = fBegin;
      mPLC.holes[0][f][1] = fEnd;
    }

    mType = starwithhole;
    this->finalize();
  }

  //------------------------------------------------------------------------
  // setTrogdor2
  // No explanation necessary
  //------------------------------------------------------------------------
  void setTrogdor2() {
    this->clear();
    const unsigned nSides = 82;
    const RealType points[164]= {5.2, 2.0, 7.0, 1.5, 7.2, 2.7, 8.0, 1.2,
				 8.5, 1.2, 8.5, 3.0, 9.5, 1.6, 10.4, 2.0,
				 9.8, 3.8, 11.4, 3.3, 11.9, 5.6, 11.0, 7.0,
				 9.8, 7.8, 7.0, 9.0, 5.6, 10.5, 4.9, 9.6,
				 3.4, 8.8, 4.0, 8.2, 4.1, 6.5, 5.0, 6.6,
				 5.2, 6.1, 5.0, 6.0, 5.1, 5.5, 5.4, 5.5,
				 5.3, 5.0, 4.0, 5.0, 3.0, 6.0, 1.5, 7.0,
				 1.1, 8.0, 0.7, 8.7, 1.3, 9.6, 0.8, 10.2,
				 0.9, 11.7, 2.0, 12.5, 3.2, 12.3, 3.8, 13.0,
				 4.6, 13.1, 5.2, 14.0, 4.0, 14.3, 3.0, 14.0,
				 2.8, 15.0, 1.9, 15.3, 1.3, 16.2, 0.5, 16.6,
				 2.2, 17.5, 3.6, 15.9, 5.5, 14.3, 6.5, 14.8,
				 7.5, 17.2, 7.9, 19.0, 9.8, 18.5, 9.7, 17.8,
				 9.8, 17.0, 9.2, 16.8, 8.8, 16.0, 8.0, 15.9,
				 7.2, 15.1, 9.6, 15.0, 12.2, 14.0, 14.5, 14.0,
				 14.9, 14.8, 15.5, 14.8, 16.0, 14.2, 15.8, 12.7,
				 12.5, 12.7, 13.1, 11.8, 15.5, 11.0, 15.4, 10.5,
				 13.5, 10.4, 12.0, 11.3, 9.9, 12.1, 8.8, 12.3,
				 8.7, 12.0, 9.0, 11.0, 11.0, 9.2, 14.0, 7.0,
				 14.8, 5.0, 14.0, 3.0, 12.0, 1.1, 10.2, 0.5,
				 8.7, 0.4, 7.2, 0.8};

    mPLC.facets.resize( nSides, vector<int>(2) );
    for (unsigned i = 0; i != nSides; ++i) {
      unsigned j = nSides - i - 1;
      mPLCpoints.push_back(points[2*j  ]);
      mPLCpoints.push_back(points[2*j+1]);
      mPLC.facets[i][0] = i;
      mPLC.facets[i][1] = (i+1) % nSides;
    }
    POLY_CHECK(mPLCpoints.size()/2 == nSides);

//     const unsigned nHolePoints = 3;
//     const RealType holePoints[6] = {15.0, 14.0, 15.3, 14.4, 15.6, 14.0};
    
//     mPLC.holes[0].resize(nHolePoints);
//     for (unsigned i = 0; i != nHolePoints; ++i) {
//       mPLCpoints.push_back(holePoints[2*i  ]);
//       mPLCpoints.push_back(holePoints[2*i+1]);
//       unsigned iBegin = mPLCpoints.size()/2 - nHolePoints + i;
//       unsigned iEnd   = mPLCpoints.size()/2 - nHolePoints + ((i + 1) % nHolePoints);
//       mPLC.holes[0][i].resize(2);
//       mPLC.holes[0][i][0] = iBegin;
//       mPLC.holes[0][i][1] = iEnd;
//     }
    
    mType = trogdor2;
    this->finalize();
  }

  //------------------------------------------------------------------------
  //-------------------- ADDITIONAL HELPER FUNCTIONS -----------------------
  //------------------------------------------------------------------------
  
  //------------------------------------------------------------------------
  // testInside
  // Tests if a given point (x,y) lies in the interior of the domain.
  // If holes are defined in mPLC, then each hole is tested separately.
  // Returns true if the point is inside the facets but outside the holes.
  //------------------------------------------------------------------------
  bool testInside( RealType* pos ) {
    Point<2, RealType> p(pos[0], pos[1]);
    return mQPLC.within(p);
    // POLY_CHECK( mPLCpoints.size() > 0 );
    // const unsigned nSides = mPLC.facets.size();
    // bool isInside = inside(x,y,nSides,offset);
    
    // for (unsigned hIt = 0; hIt != mPLC.holes.size(); ++hIt ) {
    //   isInside ^= inside(x,y,mPLC.holes[hIt].size(),offset);
    // }
    // return isInside;
  }
  
  
  //------------------------------------------------------------------------
  // inside
  // Tests if (x,y) is inside the nSide-sided polygon defined by the ordered
  // set of points in mPLCpoints starting at index 'offset'
  //------------------------------------------------------------------------
  bool inside( const RealType x, const RealType y, 
	       const unsigned nSides, unsigned& offset ) {
    unsigned j = nSides - 1;
    bool isInside = false;
    for (unsigned i = 0; i < nSides; ++i ) {
      unsigned ix = 2*(i+offset),   iy = 2*(i+offset)+1;
      unsigned jx = 2*(j+offset),   jy = 2*(j+offset)+1;
      if( ((mPLCpoints[iy] <  y  &&  mPLCpoints[jy] >= y)  ||
	   (mPLCpoints[jy] <  y  &&  mPLCpoints[iy] >= y)) &&
	  (mPLCpoints[ix] <= x  ||  mPLCpoints[jx] <= x) )
	{
	  isInside ^= ( mPLCpoints[ix] + ( y         - mPLCpoints[iy] ) /
			( mPLCpoints[jy] - mPLCpoints[iy] ) *
			( mPLCpoints[jx] - mPLCpoints[ix] ) < x );
	}
      j = i;
    }
    offset += nSides;
    return isInside;
  }

  //------------------------------------------------------------------------
  // getBoundingCircle
  // Get the maximal L-2 norm of the boundary generator set about center pt
  // NOTE: method is general to 2D or 3D. Dimension is set to 2 here
  //------------------------------------------------------------------------
  void getBoundingRadius(RealType& radius) {
    POLY_CHECK( mCenter != 0 );
    radius = 0;
    for (unsigned i = 0; i < mPLCpoints.size()/2; ++i ){
      RealType distance = 0;
      for (unsigned n = 0; n < 2; ++n ){
	distance += (mPLCpoints[2*i+n] - mCenter[n]) *
	  (mPLCpoints[2*i+n] - mCenter[n]);
      }
      radius = max( radius, sqrt( distance ) );
    }
  }

  //------------------------------------------------------------------------
  // BoostMyBoundary
  // Store the boundary info as a Boost.Geometry polygon
  //------------------------------------------------------------------------
  void boostMyBoundary() {
    boost::geometry::clear(mBGboundary);
    mBGboundary = makeBGPolygon<RealType>(mPLC, mPLCpoints);
  }

  //------------------------------------------------------------------------
  // getPointInside
  // Computes a random point 
  //------------------------------------------------------------------------  
  void getPointInside(RealType* point) {
    using IntType = typename HashKey<2>::IntType;
    using IntPoint = Point<2, IntType>;
    bool inside = false;
    IntPoint p;
    while( !inside ){
      p.x = static_cast<IntType>(random01())*mQ.maxCoord.x;
      p.y = static_cast<IntType>(random01())*mQ.maxCoord.y;
      inside = mQPLC.within(p);
    }
    Point<2, RealType> pd = mQ.dequantize(p);
    point[0] = pd.x;
    point[1] = pd.y;
  }
};
}
#endif
