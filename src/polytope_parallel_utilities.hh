#ifndef POLYTOPE_PARALLEL_UTILITIES_HH
#define POLYTOPE_PARALLEL_UTILITIES_HH
//------------------------------------------------------------------------------
// A semi-random collection of stuff that is helpful for polytope in parallel
// runs.  We just assume MPI is available here, so don't include this header
// unless that is so!
//------------------------------------------------------------------------------
#include <vector>
#include <limits>

#include <mpi.h>

#include "KeyTraits.hh"

namespace polytope {

//-----------------------------------------------------------------------------
// Traits to handle mapping RealType -> MPI data type.
//-----------------------------------------------------------------------------
template<typename RealType> struct DataTypeTraits;

template<> struct DataTypeTraits<unsigned> { 
  static MPI_Datatype MpiDataType() { return MPI_UNSIGNED; }
};

template<> struct DataTypeTraits<float> { 
  static MPI_Datatype MpiDataType() { return MPI_FLOAT; }
};

template<> struct DataTypeTraits<double> { 
  static MPI_Datatype MpiDataType() { return MPI_DOUBLE; }
};

//------------------------------------------------------------------------------
// A convenient interface to MPI_Allreduce.
//------------------------------------------------------------------------------
template<typename Value>
inline
Value
allReduce(const Value& value, 
          const MPI_Op op, 
          const MPI_Comm comm) {
  Value tmp = value;
  Value result;
  MPI_Allreduce(&tmp, &result, 1, DataTypeTraits<Value>::MpiDataType(), op, comm);
  return result;
}

}

#endif
