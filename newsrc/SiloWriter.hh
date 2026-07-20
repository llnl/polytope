#ifndef __Polytope_SiloWriter__
#define __Polytope_SiloWriter__

#include <string>
#include <float.h>
#include <map>
#include <vector>

namespace polytope {

template<int Dimension, typename RealType> class Tessellation;

//! \class SiloWriter
//! This class provides a static interface for writing Silo files 
//! containing tessellations made by polytope.
template <int Dimension, typename RealType>
class SiloWriter
{
  // No general recipe
};

//! Partial specialization for 2D tessellations.
template <typename RealType>
class SiloWriter<2, RealType>
{
  public:

  //! Write an arbitrary polygonal mesh, an associated set of 
  //! (node, edge, face, cell)-centered fields, and a corresponding set of 
  //! tags, to a SILO file in the given directory.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void write(const Tessellation<2, RealType>& mesh, 
                    const std::map<std::string, RealType*>& nodeFields,
                    const std::map<std::string, std::vector<int>*>& nodeTags,
                    const std::map<std::string, RealType*>& edgeFields,
                    const std::map<std::string, std::vector<int>*>& edgeTags,
                    const std::map<std::string, RealType*>& faceFields,
                    const std::map<std::string, std::vector<int>*>& faceTags,
                    const std::map<std::string, RealType*>& cellFields,
                    const std::map<std::string, std::vector<int>*>& cellTags,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int cycle,
                    RealType time,
                    int numFiles = -1,
                    int mpiTag = 0);

  //! Write an arbitrary polygonal mesh and an associated set of 
  //! (node, edge, face, cell)-centered fields to a SILO file in the given directory.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void write(const Tessellation<2, RealType>& mesh, 
                    const std::map<std::string, RealType*>& nodeFields,
                    const std::map<std::string, RealType*>& edgeFields,
                    const std::map<std::string, RealType*>& faceFields,
                    const std::map<std::string, RealType*>& cellFields,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int cycle,
                    RealType time,
                    int numFiles = -1,
                    int mpiTag = 0)
  {
    // Just call the general function with no tags.
    std::map<std::string, std::vector<int>*> nodeTags, edgeTags, faceTags, cellTags;
    write(mesh, nodeFields, nodeTags, edgeFields, edgeTags, faceFields, faceTags, 
          cellFields, cellTags, filePrefix, directory, cycle, time, numFiles, mpiTag);
  }

  //! Write an arbitrary polygonal mesh and an associated set of 
  //! (node, edge, face, cell)-centered fields to a SILO file. This version generates a 
  //! directory name automatically. For parallel runs, the directory 
  //! name is filePrefix-nproc. For serial runs, the directory is 
  //! the current working directory.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void write(const Tessellation<2, RealType>& mesh, 
                    const std::map<std::string, RealType*>& nodeFields,
                    const std::map<std::string, RealType*>& edgeFields,
                    const std::map<std::string, RealType*>& faceFields,
                    const std::map<std::string, RealType*>& cellFields,
                    const std::string& filePrefix,
                    int cycle,
                    RealType time,
                    int numFiles = -1,
                    int mpiTag = 0)
  {
    write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, "", cycle, time, numFiles, mpiTag);
  }

  //! This version of write omits the cycle and time arguments.
  static void write(const Tessellation<2, RealType>& mesh, 
                    const std::map<std::string, RealType*>& nodeFields,
                    const std::map<std::string, RealType*>& edgeFields,
                    const std::map<std::string, RealType*>& faceFields,
                    const std::map<std::string, RealType*>& cellFields,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int numFiles = -1,
                    int mpiTag = 0)
  {
    write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, directory, -1, -FLT_MAX,
          numFiles, mpiTag);
  }

  //! This version of write omits the cycle and time arguments and 
  //! automatically generates the directory name from the file prefix.
  static void write(const Tessellation<2, RealType>& mesh, 
                    const std::map<std::string, RealType*>& nodeFields,
                    const std::map<std::string, RealType*>& edgeFields,
                    const std::map<std::string, RealType*>& faceFields,
                    const std::map<std::string, RealType*>& cellFields,
                    const std::string& filePrefix,
                    int numFiles = -1,
                    int mpiTag = 0)
  {
    write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, "", -1, -FLT_MAX,
          numFiles, mpiTag);
  }

  static void writePoint(const std::vector<RealType>& pointCoords,
                         const std::string& filename,
                         const std::string& name,
                         const std::string& dirname,
                         int cycle,
                         RealType time);

};

//! Partial specialization for 3D tessellations.
template <typename RealType>
class SiloWriter<3, RealType>
{
  public:

  //! Write an arbitrary polygonal mesh, an associated set of 
  //! (node, edge, face, cell)-centered fields, and a corresponding set of 
  //! tags, to a SILO file in the given directory.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void write(const Tessellation<3, RealType>& mesh, 
                    const std::map<std::string, RealType*>& nodeFields,
                    const std::map<std::string, std::vector<int>*>& nodeTags,
                    const std::map<std::string, RealType*>& edgeFields,
                    const std::map<std::string, std::vector<int>*>& edgeTags,
                    const std::map<std::string, RealType*>& faceFields,
                    const std::map<std::string, std::vector<int>*>& faceTags,
                    const std::map<std::string, RealType*>& cellFields,
                    const std::map<std::string, std::vector<int>*>& cellTags,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int cycle,
                    RealType time,
                    int numFiles = -1,
                    int mpiTag = 0);

  //! Write an arbitrary polyhedral mesh and an associated set of 
  //! (node, edge, face, cell)-centered fields to a SILO file in the given directory.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void write(const Tessellation<3, RealType>& mesh, 
                    const std::map<std::string, RealType*>& nodeFields,
                    const std::map<std::string, RealType*>& edgeFields,
                    const std::map<std::string, RealType*>& faceFields,
                    const std::map<std::string, RealType*>& cellFields,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int cycle,
                    RealType time,
                    int numFiles = -1,
                    int mpiTag = 0)
  {
    // Just call the general function with no tags.
    std::map<std::string, std::vector<int>*> nodeTags, edgeTags, faceTags, cellTags;
    write(mesh, nodeFields, nodeTags, edgeFields, edgeTags, faceFields, faceTags, 
          cellFields, cellTags, filePrefix, directory, cycle, time, numFiles, mpiTag);
  }

  //! Write an arbitrary polyhedral mesh and an associated set of 
  //! (node, edge, face, cell)-centered fields to a SILO file. This version generates a 
  //! directory name automatically. For parallel runs, the directory 
  //! name is filePrefix-nproc. For serial runs, the directory is 
  //! the current working directory.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void write(const Tessellation<3, RealType>& mesh, 
                    const std::map<std::string, RealType*>& nodeFields,
                    const std::map<std::string, RealType*>& edgeFields,
                    const std::map<std::string, RealType*>& faceFields,
                    const std::map<std::string, RealType*>& cellFields,
                    const std::string& filePrefix,
                    int cycle,
                    RealType time,
                    int numFiles = -1,
                    int mpiTag = 0)
  {
    write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, "", cycle, time, numFiles, mpiTag);
  }

  //! This version of write omits the cycle and time arguments.
  static void write(const Tessellation<3, RealType>& mesh, 
                    const std::map<std::string, RealType*>& nodeFields,
                    const std::map<std::string, RealType*>& edgeFields,
                    const std::map<std::string, RealType*>& faceFields,
                    const std::map<std::string, RealType*>& cellFields,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int numFiles = -1,
                    int mpiTag = 0)
  {
    write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, directory, -1, -FLT_MAX,
          numFiles, mpiTag);
  }

  //! This version of write omits the cycle and time arguments and 
  //! automatically generates the directory name from the file prefix.
  static void write(const Tessellation<3, RealType>& mesh, 
                    const std::map<std::string, RealType*>& nodeFields,
                    const std::map<std::string, RealType*>& edgeFields,
                    const std::map<std::string, RealType*>& faceFields,
                    const std::map<std::string, RealType*>& cellFields,
                    const std::string& filePrefix,
                    int numFiles = -1,
                    int mpiTag = 0)
  {
    write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, "", -1, -FLT_MAX,
          numFiles, mpiTag);
  }

};

} // end namespace

#endif
