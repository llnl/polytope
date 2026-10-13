#ifndef __Polytope_SiloReader__
#define __Polytope_SiloReader__

#include <string>
#include <float.h>
#include <map>
#include <utility>
#include <vector>
#include "SiloUtils.hh"

namespace polytope {

namespace Silo {

// Helper function for finding available cycles.
std::vector<int> findAvailableCycles(const std::string& prefix,
                                     const std::string& directory);

}

//! \class SiloReader
//! This class provides a static interface for reading Silo files 
//! containing tessellations made by polytope.
template <int Dimension, typename TessType>
class SiloReader {
  // No general recipe
};

//! Partial specialization for 2D tessellations.
template <typename TessType>
class SiloReader<2, TessType> {
public:
  using FieldMap = std::map<std::string, std::vector<double>>;
  using FieldTypeMap = std::map<FieldCentering, FieldMap>;

  //! Returns a list of cycle numbers for Silo files dumped by a SiloWriter
  //! with the given prefix, in the given directory. If the directory is 
  //! omitted, its name is generated automatically from the prefix.
  static std::vector<int> availableCycles(const std::string& filePrefix,
                                          const std::string& directory = "") {
    return Silo::findAvailableCycles(filePrefix, directory);
  }

  //! Read an arbitrary polygonal mesh and associated field/t maps from a
  //! SILO file. Field maps are keyed by Polytope field centering.
  //! \param fields A map that will store arrays of field data read in from 
  //!               the file. If \a fields contains keys, only those fields
  //!               with those keys will be read from the file, and an error 
  //!               will occur if any of the keys are not found. If \a fields 
  //!               is empty, all data will be read in from the file.
  static void read(TessType& mesh,
                   FieldTypeMap& fields,
                   const std::string& masterFilename);

};

//! Partial specialization for 3D tessellations.
template <typename TessType>
class SiloReader<3, TessType> {
public:
  using FieldMap = std::map<std::string, std::vector<double>>;
  using FieldTypeMap = std::map<FieldCentering, FieldMap>;

  //! Returns a list of cycle numbers for Silo files dumped by a SiloWriter
  //! with the given prefix, in the given directory. If the directory is 
  //! omitted, its name is generated automatically from the prefix.
  static std::vector<int> availableCycles(const std::string& filePrefix,
                                          const std::string& directory = "") {
    return Silo::findAvailableCycles(filePrefix, directory);
  }

  //! Read an arbitrary polyhedral mesh and an associated set of 
  //! fields from a SILO file. Field maps are keyed by Polytope
  //! field centering.
  //! \param fields A map that will store arrays of field data read in from 
  //!               the file. If \a fields contains keys, only those fields
  //!               with those keys will be read from the file, and an error 
  //!               will occur if any of the keys are not found. If \a fields 
  //!               is empty, all data will be read in from the file.
  static void read(TessType& mesh,
                   FieldTypeMap& fields,
                   const std::string& masterFilename);

};

} // end namespace

#endif
