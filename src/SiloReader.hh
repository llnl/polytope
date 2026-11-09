#ifndef __Polytope_SiloReader__
#define __Polytope_SiloReader__

#include <string>
#include <vector>
#include "SiloUtils.hh"

namespace polytope {

namespace Silo {

// Helper function for finding available cycles.
std::vector<int> findAvailableCycles(const std::string& prefix,
                                     const std::string& directory);

}

//! \class SiloReader
//! This class provides a interface for reading Silo files
//! containing tessellations made by polytope.
template <int Dimension, typename TessType>
class SiloReader {
  // No general recipe
};

//! Partial specialization for 2D tessellations.
template <typename TessType>
class SiloReader<2, TessType> {
public:
  //! Returns a list of cycle numbers for Silo files dumped by a SiloWriter
  //! with the given prefix, in the given directory. If the directory is
  //! omitted, its name is generated automatically from the prefix.
  std::vector<int> availableCycles(const std::string& filePrefix,
                                   const std::string& directory = "") {
    return Silo::findAvailableCycles(filePrefix, directory);
  }

  //! Read an arbitrary polygonal mesh and all its fields from a Silo file.
  //! The fields are stored in mesh.fields.
  void read(TessType& mesh,
            const std::string& masterFilename);
};

//! Partial specialization for 3D tessellations.
template <typename TessType>
class SiloReader<3, TessType> {
public:
  //! Returns a list of cycle numbers for Silo files dumped by a SiloWriter
  //! with the given prefix, in the given directory. If the directory is
  //! omitted, its name is generated automatically from the prefix.
  std::vector<int> availableCycles(const std::string& filePrefix,
                                   const std::string& directory = "") {
    return Silo::findAvailableCycles(filePrefix, directory);
  }

  //! Read an arbitrary polyhedral mesh and all its fields from a Silo file.
  //! The fields are stored in mesh.fields.
  void read(TessType& mesh,
            const std::string& masterFilename);
};

} // end namespace

#endif
