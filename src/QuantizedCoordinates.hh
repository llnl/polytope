#ifndef POLYTOPE_QUANTIZEDCOORDINATES_HH
#define POLYTOPE_QUANTIZEDCOORDINATES_HH

#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include "ReducedPLC.hh"
#include "Point.hh"
#include "DimensionTraits.hh"
#include "polytope_geometric_utilities.hh"

namespace polytope {


//! \class QuantizedCoordinates - A Piecewise Linear Complex in 3D, or a 
//! Planar Straight Line Graph (PSLG) in 2D.
//! The quantized coordinates is a reduced PLC which
//!  a) stores lower and upper bounds for point quantization
//!  b) stores the center and radius of an inscribed bounding sphere
//!  c) has routines for passing points to/from floating point from/to integers
template<int Dimension, typename RealType>
class QuantizedCoordinates: public ReducedPLC<Dimension, RealType> {
public:
  //-------------- Public types and member variables ----------------- //

  typedef typename DimensionTraits<Dimension, RealType>::CoordHash CoordHash;
  typedef typename DimensionTraits<Dimension, RealType>::IntPoint  IntPoint;
  typedef typename DimensionTraits<Dimension, RealType>::RealPoint RealPoint;

  //! The bounding box size
  RealPoint low_inner, high_inner;
  RealPoint low_outer, high_outer;

  //! The degeneracy spacing
  RealType delta;

private:
  //-------------- Private types and member variables ----------------- //

  //! Infinite bounding sphere
  RealType mRinf_inner, mRinf_outer;
  RealPoint mCenter;

  //! Parameters for controlling degeneracy spacing
  RealType mDegeneracy;
  CoordHash mCoordMax;

  //! Flag for checking if coordinates have been modified
  bool mCoordinatesModified;

   //! Flag for checking if two bounding boxes are being used
  bool mUsingTwoBoundingBoxes;

  //! Verbosity
  bool mVerbose;

public:
  //-------------------- Public interface ---------------------- //

  //------------------------------------------------------------------------
  //! Constructors
  //------------------------------------------------------------------------

  //! Default constructor
  QuantizedCoordinates(const bool verbose = false):
    ReducedPLC<Dimension, RealType>(),
    low_inner(), 
    high_inner(), 
    low_outer(), 
    high_outer(), 
    mCenter(),
    mCoordinatesModified(false),
    mUsingTwoBoundingBoxes(false),
    mVerbose(verbose) { this->clear(); }
   
  //! Constuctor + initialize
  QuantizedCoordinates(const std::vector<RealType>& points,
                       const bool verbose = false):
    ReducedPLC<Dimension, RealType>(),
    low_inner(), 
    high_inner(),
    low_outer(), 
    high_outer(),
    mCenter(),
    mCoordinatesModified(false),
    mUsingTwoBoundingBoxes(false),
    mVerbose(verbose)
  {
    initialize(points);
    POLY_ASSERT(mRinf_inner   != 0.0);
    POLY_ASSERT(delta  != 0.0);
    POLY_ASSERT(not mCoordinatesModified);
    POLY_ASSERT(not mUsingTwoBoundingBoxes);
    storeBoundsAsPLCData();
  }

  //------------------------------------------------------------------------
  //! Clear and empty operations
  //------------------------------------------------------------------------
  void clear() {
    low_inner = RealPoint();
    high_inner = RealPoint();
    low_outer = RealPoint();
    high_outer = RealPoint();
    mCenter = RealPoint();
    mRinf_inner = 0; 
    mRinf_outer = 0; 
    delta = 0;
    mCoordinatesModified=false;
    mUsingTwoBoundingBoxes=false;
    this->points.clear();
    this->facets.clear();
    this->holes.clear();
  }

  //------------------------------------------------------------------------
  //! Degeneracy accessors
  //------------------------------------------------------------------------
  RealType degeneracy() { return mDegeneracy; }
  void degeneracy(const RealType degeneracy) {
    mDegeneracy = degeneracy;
    mCoordinatesModified = true;
  }

  //------------------------------------------------------------------------
  //! CoordMax accessor
  //------------------------------------------------------------------------
  CoordHash coordMax() const { return mCoordMax; }

  //------------------------------------------------------------------------
  //! Infinite Radius accessor
  //------------------------------------------------------------------------
  RealType infiniteRadius() const { return std::max(mRinf_inner, mRinf_outer); }

  //------------------------------------------------------------------------
  //! Center point accessor
  //------------------------------------------------------------------------
  RealPoint center() const { return mCenter; }

  //------------------------------------------------------------------------
  //! Initialize the coordinate system
  //------------------------------------------------------------------------
  void initialize(const std::vector<RealType>& points,
                  const RealType degeneracy) {
    this->clear();

    // Initialize degeneracy and coordMax
    if (not mCoordinatesModified) {
      mDegeneracy = degeneracy;
      delta = degeneracy;
    }
    
    // compute bounds for the box
    geometry::computeBoundingBox<Dimension,RealType>(points, 
                                                     true, 
                                                     &low_inner.x, 
                                                     &high_inner.x);

    // Compute the coordinate center
    computeCenterFromBounds(&low_inner.x, &high_inner.x, &mCenter.x);

    // 

    // The infinite-radius bounding sphere
    RealPoint box;
    mRinf_inner = 0.0;
    for (unsigned j = 0; j != Dimension; ++j) {
      POLY_ASSERT(low_inner[j] <= high_inner[j]);
      mCenter[j]  = 0.5*(high_inner[j] + low_inner[j]);
      box[j]      =      high_inner[j] - low_inner[j];
      mRinf_inner = std::max(mRinf_inner, 0.75*box[j]);
      POLY_ASSERT(mCenter[j] - mRinf_inner < low_inner [j] and
                  mCenter[j] + mRinf_inner > high_inner[j]);
    }

    // Resize bounds to circumscribe sphere
    computeBoundsFromRadius(mRinf_inner, &mCenter.x, &low_inner.x, &high_inner.x);

    // Figure out the maximum integer coordinate needed
    setMaximumCoordinate(2.0*mRinf_inner, mDegeneracy);
    
    mCoordinatesModified = false;
    storeBoundsAsPLCData();
    low_outer = low_inner;
    high_outer = high_inner;
    mRinf_outer = mRinf_inner;
  }

  //------------------------------------------------------------------------
  //! Expand bounding box based on new points
  //------------------------------------------------------------------------
  void expand(const std::vector<RealType>& points,
              const bool globalReduce,
              const bool twoBoundingBoxes = false) {
    RealType plow[Dimension], phigh[Dimension];
    geometry::computeBoundingBox<Dimension,RealType>(points, 
						     globalReduce, 
						     plow, 
                                                     phigh);
    expand(plow, phigh, twoBoundingBoxes);
  }

  //------------------------------------------------------------------------
  //! Expand bounding box based on new bounds
  //------------------------------------------------------------------------
  void expand(const RealType* plow, 
              const RealType* phigh,
              const bool twoBoundingBoxes = false) {
    // Pre-conditions
    POLY_ASSERT(not this->points.empty() and 
                not this->facets.empty());
    POLY_ASSERT(mRinf_outer != 0.0);
    POLY_BEGIN_CONTRACT_SCOPE;
    for (unsigned j = 0; j < Dimension; ++j) POLY_ASSERT(plow[j] <= phigh[j]);
    POLY_END_CONTRACT_SCOPE;
    
    // Resize the box
    RealType boxHalf = 0.0;
    RealPoint low_new, high_new;
    RealType rinf_new;
    for (unsigned j = 0; j != Dimension; ++j) {
      POLY_ASSERT(plow [j] <= phigh[j]);
      POLY_ASSERT(mCenter[j] <= std::max(phigh[j], high_outer[j]));
      POLY_ASSERT(mCenter[j] >= std::min(plow [j], low_outer [j]));
      boxHalf = std::max(boxHalf, std::max(phigh[j], high_outer[j]) - mCenter[j]);
      boxHalf = std::max(boxHalf, mCenter[j] - std::min(plow [j], low_outer [j]));
      // low_new [j] = std::min(low_outer [j], plow [j]);
      // high_new[j] = std::max(high_outer[j], phigh[j]);
    }
    
    // Circumscribe the new box dimensions with the infinite sphere
    rinf_new = std::max(mRinf_outer, 3.0*boxHalf);
    
    // Inflate box dimensions to inscribe the new infinite sphere 
    computeBoundsFromRadius(rinf_new, &mCenter.x, &low_new.x, &high_new.x);

    // All-reduce to ensure bounds are exact across processors.
#ifdef POLYTOPE_ENABLE_MPI
    for (unsigned j = 0; j != Dimension; ++j) {
      low_new [j] = allReduce(low_new [j], MPI_MIN, MPI_COMM_WORLD);
      high_new[j] = allReduce(high_new[j], MPI_MAX, MPI_COMM_WORLD);
    }
    rinf_new = 0.5*(high_new[0] - low_new[0]);
#endif

    // We're using two bounding boxes. Store the new bounds.
    if (twoBoundingBoxes) {
      low_outer = low_new;
      high_outer = high_new;
      mRinf_outer = rinf_new;
      mUsingTwoBoundingBoxes = true;
    }
    // Only one bounding box. The new bounds are inner. 
    // Set outer equal to inner and recompute the maximum coordinate
    else {
      low_inner = low_new;
      high_inner = high_new;
      mRinf_inner = rinf_new;
      low_outer  = low_inner;
      high_outer = high_inner;
      mRinf_outer = mRinf_inner;
      mUsingTwoBoundingBoxes = false;

      setMaximumCoordinate(2.0*mRinf_inner, mDegeneracy);
    }
  
    mCoordinatesModified = false;

    // Post conditions
    POLY_ASSERT(mRinf_outer > 0);
  }

  //------------------------------------------------------------------------
  //! Print out the bounding box info
  //------------------------------------------------------------------------
  friend std::ostream& operator<<(std::ostream& os, const QuantizedCoordinates& coords) {
     os << dynamic_cast<const ReducedPLC<Dimension, RealType>&>(coords) << std::endl
        << "Inner Bounding Box  = ";
    for (unsigned j = 0; j != Dimension; ++j) {
      os << "(" << coords.low_inner[j] << "," << coords.high_inner[j] << ") ";
    }
    os << std::endl;
    os << "Outer Bounding Box  = ";
    for (unsigned j = 0; j != Dimension; ++j) {
      os << "(" << coords.low_outer[j] << "," << coords.high_outer[j] << ") ";
    }
    os << std::endl;
    os << "Sphere Center       = " << coords.center() << std::endl;
    os << "Inf Sphere Radius   = " << coords.infiniteRadius() << std::endl;
    os << "Int Spacing         = " << coords.delta << std::endl;
    return os;
  }

  //------------------------------------------------------------------------
  //! Check if a point is inside inner bounding box
  //------------------------------------------------------------------------
  inline
  bool isInner(const RealType* point) const {
    bool result = true;
    if (not mUsingTwoBoundingBoxes) return result;
    else {
      for (unsigned j = 0; j < Dimension; ++j) {
        result *= ((low_inner[j] <= point[j]) and (point[j] <= high_inner[j]));
      }
      return result;
    }
  }
   
  //------------------------------------------------------------------------
  //! Quantize a floating-point-precision point
  //------------------------------------------------------------------------
  inline
  IntPoint quantize(const RealType* pointIn) const {
    POLY_ASSERT(this->isInner(pointIn));
    POLY_ASSERT(not mCoordinatesModified);
    return DimensionTraits<Dimension, RealType>::constructPoint(pointIn, &low_inner.x, delta, 0);
  }
  
  //------------------------------------------------------------------------
  //! Dequantize an integer point to floating-point-precision 
  //------------------------------------------------------------------------
  inline
  RealPoint dequantize(const CoordHash* pointIn) const {
    POLY_ASSERT(not mCoordinatesModified);
    RealType p[Dimension];
    for (unsigned j = 0; j < Dimension; ++j) {
      p[j] = low_inner[j] + (RealType(pointIn[j]))*delta;
      // p[j] = low_inner[j] + (RealType(pointIn[j]) + 0.5)*delta;
    }
    POLY_ASSERT(isInner(p));
    return DimensionTraits<Dimension, RealType>::constructPoint(p);
  }

  //------------------------------------------------------------------------
  //! Project a point to the bounding circle
  //------------------------------------------------------------------------
  inline
  RealPoint projectPoint(const RealType* rayPoint,
                         const RealType* rayDirection) const {
    RealPoint result;
    bool test = geometry::rayCircleIntersection(rayPoint, rayDirection, 
                                                &mCenter.x, mRinf_outer, 
                                                delta, &result.x);
    POLY_CONTRACT_VAR(test);
    POLY_ASSERT(test);
    return result;
  }

private:
  //-------------------- Private interface ---------------------- //

  //------------------------------------------------------------------------
  //! Fill in ReducedPLC info
  //------------------------------------------------------------------------
  void storeBoundsAsPLCData() {
    std::vector<RealType> pts;
    std::vector<std::vector<int> > fcts;
    if (Dimension == 2) {
      pts.push_back(low_inner [0]);  pts.push_back(low_inner [1]);
      pts.push_back(high_inner[0]);  pts.push_back(low_inner [1]);
      pts.push_back(high_inner[0]);  pts.push_back(high_inner[1]);
      pts.push_back(low_inner [0]);  pts.push_back(high_inner[1]);

      fcts.resize(4, std::vector<int>(2));
      for (unsigned j = 0; j != 4; ++j) {
        fcts[j][0] = j;  fcts[j][1] = (j+1)%4;
      }
    } else {
      POLY_ASSERT(Dimension == 3);
      pts.push_back(low_inner [0]);  pts.push_back(low_inner [1]);  pts.push_back(low_inner [2]);
      pts.push_back(high_inner[0]);  pts.push_back(low_inner [1]);  pts.push_back(low_inner [2]);
      pts.push_back(high_inner[0]);  pts.push_back(high_inner[1]);  pts.push_back(low_inner [2]);
      pts.push_back(low_inner [0]);  pts.push_back(high_inner[1]);  pts.push_back(low_inner [2]);
      pts.push_back(low_inner [0]);  pts.push_back(low_inner [1]);  pts.push_back(high_inner[2]);
      pts.push_back(high_inner[0]);  pts.push_back(low_inner [1]);  pts.push_back(high_inner[2]);
      pts.push_back(high_inner[0]);  pts.push_back(high_inner[1]);  pts.push_back(high_inner[2]);
      pts.push_back(low_inner [0]);  pts.push_back(high_inner[1]);  pts.push_back(high_inner[2]);

      fcts.resize(6, std::vector<int>(4));
      fcts[0][0] = 0;  fcts[0][1] = 1;  fcts[0][2] = 2;  fcts[0][3] = 3;  // -z face
      fcts[0][0] = 4;  fcts[0][1] = 5;  fcts[0][2] = 6;  fcts[0][3] = 7;  // +z face
      fcts[0][0] = 0;  fcts[0][1] = 1;  fcts[0][2] = 5;  fcts[0][3] = 4;  // -y face
      fcts[0][0] = 2;  fcts[0][1] = 6;  fcts[0][2] = 7;  fcts[0][3] = 3;  // +y face
      fcts[0][0] = 0;  fcts[0][1] = 4;  fcts[0][2] = 7;  fcts[0][3] = 3;  // -x face
      fcts[0][0] = 1;  fcts[0][1] = 5;  fcts[0][2] = 6;  fcts[0][3] = 2;  // +x face
    }
    this->points = pts;
    this->facets = fcts;
  }   

  //------------------------------------------------------------------------
  //! Compute the center point from bounds
  //------------------------------------------------------------------------
  void computeCenterFromBounds(const RealType* low,
                               const RealType* high,
                               RealType* center) {
    for (int j = 0; j < Dimension; ++j) {
      POLY_ASSERT(low[j] <= high[j]);
      center[j] = 0.5*(low[j] + high[j]);
    }
  }

  //------------------------------------------------------------------------
  //! Set bounds that circumscribe a given circle
  //------------------------------------------------------------------------
  void computeRadiusFromBounds(const RealType* low,
                               const RealType* high,
                               const RealType* center,
                               const RealType resizeFactor,
                               RealType radius) {
    POLY_ASSERT(resizeFactor > 0.0);
    RealType boxHalf = 0.0;
    for (unsigned j = 0; j != Dimension; ++j) {
      POLY_ASSERT(low[j] <= high[j]);
      POLY_ASSERT(center[j] <= high[j]);
      POLY_ASSERT(center[j] >= low [j]);
      boxHalf = std::max(boxHalf, high  [j] - center[j]);
      boxHalf = std::max(boxHalf, center[j] - low   [j]);
    }
    radius = resizeFactor*boxHalf;
    POLY_ASSERT(radius > 0.0);
  }

  //------------------------------------------------------------------------
  //! Set bounds that circumscribe a given circle
  //------------------------------------------------------------------------
  void computeBoundsFromRadius(const RealType radius,
                               const RealType* center,
                               RealType* low,
                               RealType* high) {
    POLY_ASSERT(radius > 0.0);
    for (unsigned j = 0; j != Dimension; ++j) {
      low [j] = std::min(low [j], center[j] - radius);
      high[j] = std::max(high[j], center[j] + radius);
      POLY_ASSERT(low[j] <= high[j]);
    }
  }

  //------------------------------------------------------------------------
  //! Set the maximum integer coordinate.
  //! NOTE: Emits a warning if the selected degeneracy cannot be reached
  //!       by 32-bit integers. Fails if 64-bit integers will not work.
  //------------------------------------------------------------------------
  void setMaximumCoordinate(const RealType size, const RealType degeneracy) {
    CoordHash maxCoord = CoordHash(ceil(size / degeneracy));
    if (maxCoord > std::numeric_limits<int32_t>::max() and mVerbose) {
       std::cerr << "WARNING: The floating point degeneracy tolerance "
                 << "specified of " << degeneracy << " overflows a 32-bit "
                 << "integer grid. Robustness issues can result when using "
                 << "Boost.Geometry to perform intersections." << std::endl;
    }
    POLY_VERIFY2(maxCoord < std::numeric_limits<CoordHash>::max(),
                 "Quantization Error: The floating point degeneracy tolerance " <<
                 "specified of " << mDegeneracy << " overflows the integer grid. " <<
                 "Please specify a coarser tolerance or use a longer integer."); 
    mCoordMax = maxCoord;
  }
};

} //end polytope namespace

#endif
