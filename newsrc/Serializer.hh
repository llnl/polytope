//----------------------------------------------------------------------------//
// A bunch of utility methods to help with serializing/unserializing objects
// in polytope.
//----------------------------------------------------------------------------//
#ifndef __Polytope_Serializer__
#define __Polytope_Serializer__

#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <string>
#include <stdint.h>

#include "polytope_internal.hh"
#include "PLC.hh"
#include "ReducedPLC.hh"

namespace polytope {

//------------------------------------------------------------------------------
// Serialize.  Due to limitations on partial specialization for functions, we 
// define a Functor object to handle specializing for each data type.
//------------------------------------------------------------------------------
// Generic primitive types.
template<typename T>
struct Serializer {

  static void serializeImpl(const T& val, 
                            std::vector<char>& buffer) {
    const unsigned n = sizeof(T);
    const char* data = reinterpret_cast<const char*>(&val);
    std::copy(data, data + n, std::back_inserter(buffer));
  }

  static void deserializeImpl(T& val, 
                              std::vector<char>::const_iterator& bufItr, 
                              const std::vector<char>::const_iterator& endItr) {
    const unsigned n = sizeof(T);
    char* data = reinterpret_cast<char*>(&val);
    POLY_CONTRACT_VAR(endItr);
    POLY_ASSERT(bufItr + n <= endItr);
    std::copy(bufItr, bufItr + n, data);
    bufItr += n;
  }
};

// std::string.
template<>
struct Serializer<std::string> {

  static void serializeImpl(const std::string& val,
                            std::vector<char>& buffer) {
    const unsigned n = val.size();
    Serializer<unsigned>::serializeImpl(n, buffer);
    std::copy(val.begin(), val.end(), std::back_inserter(buffer));
  }

  static void deserializeImpl(std::string& val,
                              std::vector<char>::const_iterator& bufItr,
                              const std::vector<char>::const_iterator& endItr) {
    unsigned n;
    Serializer<unsigned>::deserializeImpl(n, bufItr, endItr);
    POLY_CONTRACT_VAR(endItr);
    POLY_ASSERT(bufItr + n <= endItr);
    val.assign(bufItr, bufItr + n);
    bufItr += n;
  }
};

// std::vector of known type.
template<typename T>
struct Serializer<std::vector<T> > {

  static void serializeImpl(const std::vector<T>& val, 
                            std::vector<char>& buffer) {
    const unsigned n = val.size();
    Serializer<unsigned>::serializeImpl(n, buffer);
    for (unsigned i = 0; i != n; ++i) Serializer<T>::serializeImpl(val[i], buffer);
  }

  static void deserializeImpl(std::vector<T>& val, 
                              std::vector<char>::const_iterator& bufItr, 
                              const std::vector<char>::const_iterator& endItr) {
    unsigned n, n0 = val.size();
    Serializer<unsigned>::deserializeImpl(n, bufItr, endItr);
    val.resize(n0 + n);
    for (unsigned i = 0; i != n; ++i) Serializer<T>::deserializeImpl(val[n0 + i], bufItr, endItr);
  }
};

// std::vector<std::vector> of known type.
template<typename T>
struct Serializer<std::vector<std::vector<T> > > {

  static void serializeImpl(const std::vector<std::vector<T> >& val, 
                            std::vector<char>& buffer) {
    unsigned n = val.size();
    Serializer<unsigned>::serializeImpl(n, buffer);
    for (unsigned i = 0; i != n; ++i) Serializer<std::vector<T> >::serializeImpl(val[i], buffer);
  }

  static void deserializeImpl(std::vector<std::vector<T> >& val, 
                              std::vector<char>::const_iterator& bufItr, 
                              const std::vector<char>::const_iterator& endItr) {
    unsigned n, n0 = val.size();
    Serializer<unsigned>::deserializeImpl(n, bufItr, endItr);
    val.resize(n0 + n);
    for (unsigned i = 0; i != n; ++i) Serializer<std::vector<T> >::deserializeImpl(val[n0 + i], bufItr, endItr);
  }
};

// std::vector<std::vector<std::vector>> of known type.
template<typename T>
struct Serializer<std::vector<std::vector<std::vector<T> > > > {

  static void serializeImpl(const std::vector<std::vector<std::vector<T> > >& val, 
                            std::vector<char>& buffer) {
    const unsigned n = val.size();
    Serializer<unsigned>::serializeImpl(n, buffer);
    for (unsigned i = 0; i != n; ++i) Serializer<std::vector<std::vector<T> > >::serializeImpl(val[i], buffer);
  }

  static void deserializeImpl(std::vector<std::vector<std::vector<T> > >& val, 
                              std::vector<char>::const_iterator& bufItr, 
                              const std::vector<char>::const_iterator& endItr) {
    unsigned n, n0 = val.size();
    Serializer<unsigned>::deserializeImpl(n, bufItr, endItr);
    val.resize(n0 + n);
    for (unsigned i = 0; i != n; ++i) Serializer<std::vector<std::vector<T> > >::deserializeImpl(val[n0 + i], bufItr, endItr);
  }
};

//------------------------------------------------------------------------------
// Dispatch to the functors!
//------------------------------------------------------------------------------
template<typename T>
void serialize(const T& val, 
               std::vector<char>& buffer) {
  Serializer<T>::serializeImpl(val, buffer);
}

template<typename T>
void deserialize(T& val,
                 std::vector<char>::const_iterator& bufItr,
                 const std::vector<char>::const_iterator& endItr) {
  Serializer<T>::deserializeImpl(val, bufItr, endItr);
}

//------------------------------------------------------------------------------
// Serialize a PLC.
//------------------------------------------------------------------------------
template<int Dimension>
struct Serializer<PLC<Dimension>>
{

  static void serializeImpl(const PLC<Dimension>& val,
                            std::vector<char>& buffer) {
    serialize(val.facets, buffer);
    serialize(val.holes, buffer);
  }

  static void deserializeImpl(PLC<Dimension>& val,
                              std::vector<char>::const_iterator& bufItr,
                              const std::vector<char>::const_iterator& endItr) {
    deserialize(val.facets, bufItr, endItr);
    deserialize(val.holes, bufItr, endItr);
  }
};

//------------------------------------------------------------------------------
// Serialize a ReducedPLC.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
struct Serializer<ReducedPLC<Dimension, RealType> > {

  static void serializeImpl(const ReducedPLC<Dimension, RealType>& val,
                            std::vector<char>& buffer) {
    const unsigned nf = val.facets.size();
    serialize(nf, buffer);
    for (unsigned i = 0; i != nf; ++i) serialize(val.facets[i], buffer);
    serialize(val.holes, buffer);
    serialize(val.points, buffer);
  }

  static void deserializeImpl(ReducedPLC<Dimension, RealType>& val,
                              std::vector<char>::const_iterator& bufItr,
                              const std::vector<char>::const_iterator& endItr) {
    unsigned nf;
    deserialize(nf, bufItr, endItr);
    val.facets.resize(nf);
    for (unsigned i = 0; i != nf; ++i) deserialize(val.facets[i], bufItr, endItr);
    deserialize(val.holes, bufItr, endItr);
    deserialize(val.points, bufItr, endItr);
  }
};

}

#endif
