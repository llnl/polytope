#ifndef POLYTOPE_BOUNDARY3D_HH
#define POLYTOPE_BOUNDARY3D_HH

//#include "within.hh"
#include "QuantPLC.hh"
#include "Shapes.hh"
#include "Boundary2D.hh"

using namespace std;
namespace polytope {
//------------------------------------------------------------------------
class Boundary3D {
public:
  using RealType = double;
  // -------------- Public member variables and routines ---------------- //
  
  // Piecewise linear construct to define the boundary facets + holes
  PLC<3> mPLC;
  // Vector of generators to define the boundary
  std::vector<RealType> mPLCpoints;
  // Ranges of bounding box
  RealType mCenter[3];
  // Define enum to keep track fo the type of boundary called for
  enum BoundaryType{
    square             = 0,
    circle             = 1,
    circlewithstarhole = 2,
    cylinderwithhole   = 3,
    starwithhole       = 4,
  };
  
  // Boundary type
  mutable BoundaryType mType;

  Quantizer<3> mQ;
  QuantPLC<3> mQPLC;
  
  //------------------------------------------------------------------------
  // Constructor, destructor
  //------------------------------------------------------------------------
  Boundary3D():
    mType(square){
    this->clear();
    this->setDefaultBoundary(mType);
  }
  
  ~Boundary3D() {};
    
  void clear() {
    mPLC.clear();
    mPLCpoints.clear();
  }
  
  void finalize() {
    mQ.init(mPLCpoints);
    mQPLC.init(mPLC, mQ, mPLCpoints);
  }

  //------------------------------------------------------------------------
  // setDefaultBoundary
  //------------------------------------------------------------------------
  void setDefaultBoundary(const int bType)
  {
    mCenter[0] = 0.0;
    mCenter[1] = 0.0;

    switch(bType){
    case square:
      this->setUnitSquare();
      break;
    case circle:
      this->setUnitCircle();
      break;
    case circlewithstarhole:
      this->setCircleWithStarHole();
      break;
    case cylinderwithhole:
      this->setCylinderHole();
      break;
    case starwithhole:
      this->setStarWithHole();
      break;
    }
   }

  // ------------------------------------------------------------
  // Box with one cylindrical hole
  // ------------------------------------------------------------
  void setCylinderHole() {
    this->clear();

    // bpoints contains all PLC vertices, outer boundary + hole vertices
    std::vector<RealType> bpoints = {
        // Box vertices 0..7
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 1.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        1.0, 0.0, 1.0,
        1.0, 1.0, 1.0,
        0.0, 1.0, 1.0
    };

    std::vector<std::vector<int>> boundaryFacets = {
        {0, 1, 2, 3},
        {4, 5, 6, 7},
        {0, 1, 5, 4},
        {1, 2, 6, 5},
        {2, 3, 7, 6},
        {3, 0, 4, 7}
    };

    // Add cylindrical hole vertices into bpoints
    const int N = 16;
    const RealType cx = 0.5;
    const RealType cy = 0.5;
    const RealType r = 0.15;
    const RealType z0 = 0.2;
    const RealType z1 = 0.8;

    int bottomStart = static_cast<int>(bpoints.size() / 3);
    for (int i = 0; i < N; ++i) {
      RealType a = 2.0 * M_PI * i / N;
      bpoints.push_back(cx + r * std::cos(a));
      bpoints.push_back(cy + r * std::sin(a));
      bpoints.push_back(z0);
    }

    int topStart = static_cast<int>(bpoints.size() / 3);
    for (int i = 0; i < N; ++i) {
      RealType a = 2.0 * M_PI * i / N;
      bpoints.push_back(cx + r * std::cos(a));
      bpoints.push_back(cy + r * std::sin(a));
      bpoints.push_back(z1);
    }

    std::vector<std::vector<std::vector<int>>> holeFacets;
    holeFacets.resize(1);

    // Side facets
    for (int i = 0; i < N; ++i) {
      int j = (i + 1) % N;
      holeFacets[0].push_back({
                               bottomStart + i,
                               bottomStart + j,
                               topStart + j,
                               topStart + i
        });
    }

    // Bottom cap
    {
      std::vector<int> cap;
      for (int i = 0; i < N; ++i) {
        cap.push_back(bottomStart + i);
      }
      holeFacets[0].push_back(cap);
    }

    // Top cap
    {
      std::vector<int> cap;
      for (int i = N - 1; i >= 0; --i) {
        cap.push_back(topStart + i);
      }
      holeFacets[0].push_back(cap);
    }
    
    mPLC.facets = boundaryFacets;
    mPLC.holes = holeFacets;
    mPLCpoints = std::move(bpoints);
    finalize();
  }

  //------------------------------------------------------------------------
  // setUnitSquare
  //------------------------------------------------------------------------
  void setUnitSquare() {
    this->clear();

    const RealType x1 = mCenter[0] - 0.5;
    const RealType x2 = mCenter[0] + 0.5;
    const RealType y1 = mCenter[1] - 0.5;
    const RealType y2 = mCenter[1] + 0.5;
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
    RealType z = 0.;
    for (int i = 0; i < 2; ++i) {
      for (unsigned b = 0; b < Nb; ++b) {
        RealType theta = 2.0*M_PI*b/(Nb+1);
        RealType x = mCenter[0] + cos(theta);
        RealType y = mCenter[1] + sin(theta);
        mPLCpoints.push_back(x);
        mPLCpoints.push_back(y);
        mPLCpoints.push_back(z);
      }
      z = 1.;
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
    RealType zloc = 0.;
    for (int i = 0; i < 2; ++i) {
      for (unsigned p = 0; p < nPoints; ++p ) {
        // For the pointy bits of the star
        theta = M_PI/2 - p*theta0;
        mPLCpoints.push_back(mCenter[0] + outerRadius*cos(theta));
        mPLCpoints.push_back(mCenter[1] + outerRadius*sin(theta));
        mPLCpoints.push_back(zloc);

        // For the concave bits of the star
        theta = M_PI/2 - p*theta0 - theta0/2.0;
        mPLCpoints.push_back(mCenter[0] + innerRadius*cos(theta));
        mPLCpoints.push_back(mCenter[1] + innerRadius*sin(theta));
        mPLCpoints.push_back(zloc);
      }
      zloc = 1.;
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
  //-------------------- ADDITIONAL HELPER FUNCTIONS -----------------------
  //------------------------------------------------------------------------
  
  //------------------------------------------------------------------------
  // testInside
  // Tests if a given point (x,y) lies in the interior of the domain.
  // If holes are defined in mPLC, then each hole is tested separately.
  // Returns true if the point is inside the facets but outside the holes.
  //------------------------------------------------------------------------
  bool testInside(RealType* pos) {
    Point<3, RealType> p(pos[0], pos[1], pos[2]);
    return mQPLC.within(p);
  }
  
  
  //------------------------------------------------------------------------
  // inside
  // Tests if (x,y) is inside the nSide-sided polygon defined by the ordered
  // set of points in mPLCpoints starting at index 'offset'
  //------------------------------------------------------------------------
  bool inside(const RealType x, const RealType y, 
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
  // NOTE: method is general to 3D or 3D. Dimension is set to 2 here
  //------------------------------------------------------------------------
  void getBoundingRadius(RealType& radius) {
    POLY_CHECK( mCenter != 0 );
    radius = 0;
    for (unsigned i = 0; i < mPLCpoints.size()/3; ++i ){
      RealType distance = 0;
      for (unsigned n = 0; n < 2; ++n ){
	distance += (mPLCpoints[3*i+n] - mCenter[n]) *
	  (mPLCpoints[3*i+n] - mCenter[n]);
      }
      radius = max( radius, sqrt( distance ) );
    }
  }

  //------------------------------------------------------------------------
  // getPointInside
  // Computes a random point 
  //------------------------------------------------------------------------  
  void getPointInside(RealType* point) {
    using IntType = typename HashKey<3>::IntType;
    using IntPoint = Point<3, IntType>;
    bool inside = false;
    IntPoint p;
    while( !inside ){
      p.x = static_cast<IntType>(random01())*mQ.maxCoord.x;
      p.y = static_cast<IntType>(random01())*mQ.maxCoord.y;
      inside = mQPLC.within(p);
    }
    Point<3, RealType> pd = mQ.dequantize(p);
    point[0] = pd.x;
    point[1] = pd.y;
  }
};
}
#endif
