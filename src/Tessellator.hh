#ifndef POLYTOPE_TESSELLATOR_HH
#define POLYTOPE_TESSELLATOR_HH

#include "QuantTessellation.hh"
#include "polytope_internal.hh"
#include "Tessellation.hh"

namespace polytope {

//! \class Tessellator - An abstract base class for objects that generate 
//! Voronoi and Voronoi-like tessellations for sets of points and/or 
//! geometries.
template<int Dimension, typename RealType>
class Tessellator {
public:

  using QuantizedTessellation = QuantTessellation<Dimension>;
  using Quant = Quantizer<Dimension>;

  //! Default constructor.
  Tessellator() = default;
  Tessellator(const Quant& Q) : m_Q(Q) {}

  //! Destructor.
  virtual ~Tessellator() {}

  //! Generate a Voronoi tessellation for the given set of generator points.
  //! The coordinates of these points are stored in point-major order and 
  //! the 0th component of the ith point appears in points[Dimension*i].
  //! \param points A (Dimension*numPoints) array containing point coordinates.
  //! \param mesh This will store the resulting tessellation.
  virtual void tessellate(const std::vector<RealType>& points,
                          Tessellation<Dimension, RealType>& mesh) const;

  //! Generate a Voronoi-like tessellation for the given set of generator 
  //! points and a description of the geometry in which they exist.
  //! The coordinates of these points are stored in point-major order and 
  //! the 0th component of the ith point appears in points[Dimension*i].
  //! This default implementation issues an error explaining that the 
  //! Tessellator does not support PLCs.
  //! \param points A (Dimension*numPoints) array containing point coordinates.
  //! \param PLCpoints A (Dimension*n) array containing point coordinates for the PLC.
  //! \param geometry A description of the geometry in Piecewise Linear Complex form.
  //! \param mesh This will store the resulting tessellation.
  virtual void tessellate(const std::vector<RealType>& points,
                          const std::vector<RealType>& PLCpoints,
                          const PLC<Dimension>& geometry,
                          Tessellation<Dimension, RealType>& mesh) const;

  //! Generate a Voronoi tessellation for the given set of generator points
  //! with a bounding box specified by \a low and \a high. Here, low[i]
  //! contains the ith coordinate for the "lower-left-near" corner of the 
  //! bounding box in 2D or 3D, and high[i] contains the corresponding 
  //! opposite corner. The coordinates of these points are stored in 
  //! point-major order and the 0th component of the ith point appears in 
  //! points[Dimension*i].
  //! \param points A (Dimension*numPoints) array containing point coordinates.
  //! \param low The coordinates of the "lower-left-near" bounding box corner.
  //! \param high The coordinates of the "upper-right-far" bounding box corner.
  //! \param mesh This will store the resulting tessellation.
  virtual void tessellate(const std::vector<RealType>& points,
                          RealType* low,
                          RealType* high,
                          Tessellation<Dimension, RealType>& mesh) const;

  //! Generate a Voronoi-like tessellation for the given set of generator 
  //! points and a description of the geometry in which they exist.
  //! The geometry description uses the ReducedPLC to combine vertex
  //! coordinates and facet topology into a single struct out of convenience.
  //! \param points A (Dimension*numPoints) array containing point coordinates.
  //! \param geometry A description of the geometry in Reduced Piecewise Linear Complex form.
  //! \param mesh This will store the resulting tessellation.
  virtual void tessellate(const std::vector<RealType>& points,
                          const ReducedPLC<Dimension, RealType>& geometry,
                          Tessellation<Dimension, RealType>& mesh) const;

  //! The following methods all return the same sort of tessellation as the 
  //! above versions, except these versions do not assume that the input
  //! generators are unique.  We allow degeneracies here, which implies a 
  //! given tessellation cell may correspond to more than one input generator.
  //! The returned vector is the same size as the input coordinates, and 
  //! indicates which tessellation cell goes with the corresponding generator
  //! coordinates.  In the case of unique input this array will simply be a 
  //! sequentially increasing array of integers.

  //! Unbounded case.
  virtual std::vector<unsigned>
  tessellateDegenerate(const std::vector<RealType>& points,
                       const RealType tol,
                       Tessellation<Dimension, RealType>& mesh) const;
  //! Bounded by a box.
  virtual std::vector<unsigned>
  tessellateDegenerate(const std::vector<RealType>& points,
                       RealType* low,
                       RealType* high,
                       const RealType tol,
                       Tessellation<Dimension, RealType>& mesh) const;

  //! Bounded by a PLC.
  virtual std::vector<unsigned>
  tessellateDegenerate(const std::vector<RealType>& points,
                       const std::vector<RealType>& PLCpoints,
                       const PLC<Dimension>& geometry,
                       const RealType tol,
                       Tessellation<Dimension, RealType>& mesh) const;

  //! Bounded by a PLC.
  virtual std::vector<unsigned>
  tessellateDegenerate(const std::vector<RealType>& points,
                       const ReducedPLC<Dimension, RealType>& geometry,
                       const RealType tol,
                       Tessellation<Dimension, RealType>& mesh) const;


  //! Override this method to return true if this Tessellator supports 
  //! the description of a domain boundary using a PLC (as in the second 
  //! tessellate method, above), and false if it does not. Some algorithms 
  //! for tessellation do not naturally accommodate an explicit boundary 
  //! description, and Tessellators using these algorithms should override 
  //! this method to return false. A stub method for PLC-enabled
  //! tessellation is provided for convenience.
  //! This query mechanism prevents us from descending into the taxonomic 
  //! hell associated with elaborate inheritance hierarchies.
  virtual bool handlesPLCs() const { return true; }

  //! Required for all tessellators:
  //! Compute the quantized tessellation.  This is the basic method all
  //! Tessellator implementations must provide, on which the other tessellation methods
  //! in polytope build.
  virtual void
  tessellateQuantized(QuantizedTessellation& qmesh) const = 0;

  //! Required for all tessellators:
  //! A unique name string per tessellation instance.
  virtual std::string name() const = 0;

  //! Required for all tessellators:
  //! Returns the accuracy to which this tessellator can distinguish coordinates.
  //! Should be returned appropriately for normalized coordinates, i.e., if all
  //! coordinates are in the range xi \in [0,1], what is the minimum allowed 
  //! delta in x.
  virtual RealType degeneracy() const = 0;
  virtual void degeneracy(const RealType val) const {};

private:

  Quant m_Q;
  // Disallowed.
  Tessellator(const Tessellator&);
  Tessellator& operator=(const Tessellator&);
};

}

#include "TessellatorInline.hh"

#endif
