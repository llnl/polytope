#ifndef __Polytope_SiloWriter__
#define __Polytope_SiloWriter__

#include <string>
#include <float.h>
#include <map>
#include <vector>
#include "Communicator.hh"

namespace polytope {

template<int Dimension, typename RealType> class Tessellation;

//! \class SiloWriter
//! This class provides a static interface for writing Silo files 
//! containing tessellations made by polytope.
template <int Dimension, typename TessType>
class SiloWriter {
  // No general recipe
};

//! Partial specialization for 2D tessellations.
template <typename TessType>
class SiloWriter<2, TessType> {
  public:

  //! Write an arbitrary polygonal mesh, an associated set of 
  //! (node, edge, face, cell)-centered fields, and a corresponding set of 
  //! tags, to a SILO file in the given directory.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void write(const TessType& mesh, 
                    const std::map<std::string, double*>& nodeFields,
                    const std::map<std::string, std::vector<int>*>& nodeTags,
                    const std::map<std::string, double*>& edgeFields,
                    const std::map<std::string, std::vector<int>*>& edgeTags,
                    const std::map<std::string, double*>& faceFields,
                    const std::map<std::string, std::vector<int>*>& faceTags,
                    const std::map<std::string, double*>& cellFields,
                    const std::map<std::string, std::vector<int>*>& cellTags,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int cycle,
                    double time,
                    int numFiles = -1,
                    int mpiTag = 0);

  //! Write an arbitrary polygonal mesh and an associated set of 
  //! (node, edge, face, cell)-centered fields to a SILO file in the given directory.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void write(const TessType& mesh, 
                    const std::map<std::string, double*>& nodeFields,
                    const std::map<std::string, double*>& edgeFields,
                    const std::map<std::string, double*>& faceFields,
                    const std::map<std::string, double*>& cellFields,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int cycle,
                    double time,
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
  static void write(const TessType& mesh, 
                    const std::map<std::string, double*>& nodeFields,
                    const std::map<std::string, double*>& edgeFields,
                    const std::map<std::string, double*>& faceFields,
                    const std::map<std::string, double*>& cellFields,
                    const std::string& filePrefix,
                    int cycle,
                    double time,
                    int numFiles = -1,
                    int mpiTag = 0)
  {
    write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, "", cycle, time, numFiles, mpiTag);
  }

  //! This version of write omits the cycle and time arguments.
  static void write(const TessType& mesh, 
                    const std::map<std::string, double*>& nodeFields,
                    const std::map<std::string, double*>& edgeFields,
                    const std::map<std::string, double*>& faceFields,
                    const std::map<std::string, double*>& cellFields,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int numFiles = -1,
                    int mpiTag = 0)
  {
    write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, directory, -1, -1.,
          numFiles, mpiTag);
  }

  //! This version of write omits the cycle and time arguments and 
  //! automatically generates the directory name from the file prefix.
  static void write(const TessType& mesh, 
                    const std::map<std::string, double*>& nodeFields,
                    const std::map<std::string, double*>& edgeFields,
                    const std::map<std::string, double*>& faceFields,
                    const std::map<std::string, double*>& cellFields,
                    const std::string& filePrefix,
                    int numFiles = -1,
                    int mpiTag = 0)
  {
    write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, "", -1, -1.,
          numFiles, mpiTag);
  }

  static void write(const TessType& mesh,
                    const std::string& filePrefix,
                    int numFiles = -1,
                    int mpiTag = 0) {
    std::map<std::string, double*> nodeFields, edgeFields, faceFields, cellFields;
    size_t meshSize = mesh.cells.size();
    std::vector<double> index(meshSize);
    std::vector<double> genx (meshSize);
    std::vector<double> geny (meshSize);
    for (int i = 0; i < meshSize; ++i) {
      index[i] = double(i);
      genx[i] = static_cast<double>(mesh.points[i].x);
      geny[i] = static_cast<double>(mesh.points[i].y);
    }
    cellFields["cell_index"] = &index[0];
    cellFields["gen_x"     ] = &genx[0];
    cellFields["gen_y"     ] = &geny[0];
#ifdef POLYTOPE_ENABLE_MPI
    int rank = Communicator::getRank();
    std::vector<double> rankField(meshSize, static_cast<double>(rank));
    cellFields["rank"      ] = &rankField[0];
#endif
    write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, "", -1, -1.,
          numFiles, mpiTag);
  }

};

//! Partial specialization for 3D tessellations.
template <typename TessType>
class SiloWriter<3, TessType> {
public:

  //! Write an arbitrary polygonal mesh, an associated set of 
  //! (node, edge, face, cell)-centered fields, and a corresponding set of 
  //! tags, to a SILO file in the given directory.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void write(const TessType& mesh, 
                    const std::map<std::string, double*>& nodeFields,
                    const std::map<std::string, std::vector<int>*>& nodeTags,
                    const std::map<std::string, double*>& edgeFields,
                    const std::map<std::string, std::vector<int>*>& edgeTags,
                    const std::map<std::string, double*>& faceFields,
                    const std::map<std::string, std::vector<int>*>& faceTags,
                    const std::map<std::string, double*>& cellFields,
                    const std::map<std::string, std::vector<int>*>& cellTags,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int cycle,
                    double time,
                    int numFiles = -1,
                    int mpiTag = 0);

  //! Write an arbitrary polyhedral mesh and an associated set of 
  //! (node, edge, face, cell)-centered fields to a SILO file in the given directory.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void write(const TessType& mesh, 
                    const std::map<std::string, double*>& nodeFields,
                    const std::map<std::string, double*>& edgeFields,
                    const std::map<std::string, double*>& faceFields,
                    const std::map<std::string, double*>& cellFields,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int cycle,
                    double time,
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
  static void write(const TessType& mesh, 
                    const std::map<std::string, double*>& nodeFields,
                    const std::map<std::string, double*>& edgeFields,
                    const std::map<std::string, double*>& faceFields,
                    const std::map<std::string, double*>& cellFields,
                    const std::string& filePrefix,
                    int cycle,
                    double time,
                    int numFiles = -1,
                    int mpiTag = 0)
  {
    write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, "", cycle, time, numFiles, mpiTag);
  }

  //! This version of write omits the cycle and time arguments.
  static void write(const TessType& mesh, 
                    const std::map<std::string, double*>& nodeFields,
                    const std::map<std::string, double*>& edgeFields,
                    const std::map<std::string, double*>& faceFields,
                    const std::map<std::string, double*>& cellFields,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int numFiles = -1,
                    int mpiTag = 0)
  {
    write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, directory, -1, -1.,
          numFiles, mpiTag);
  }

  //! This version of write omits the cycle and time arguments and 
  //! automatically generates the directory name from the file prefix.
  static void write(const TessType& mesh, 
                    const std::map<std::string, double*>& nodeFields,
                    const std::map<std::string, double*>& edgeFields,
                    const std::map<std::string, double*>& faceFields,
                    const std::map<std::string, double*>& cellFields,
                    const std::string& filePrefix,
                    int numFiles = -1,
                    int mpiTag = 0)
  {
    write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, "", -1, -1.,
          numFiles, mpiTag);
  }

  static void write(const TessType& mesh,
                    const std::string& filePrefix,
                    int numFiles = -1,
                    int mpiTag = 0) {
    std::map<std::string, double*> nodeFields, edgeFields, faceFields, cellFields;
    write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, "", -1, -1.,
          numFiles, mpiTag);
  }

};

} // end namespace

#endif
