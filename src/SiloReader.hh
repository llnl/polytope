#ifndef POLYTOPE_SILO_READER_HH
#define POLYTOPE_SILO_READER_HH

#ifdef POLYTOPE_ENABLE_MPI
#include <mpi.h>
#define MMPI_Comm MPI_Comm
#define MMPI_COMM_WORLD MPI_COMM_WORLD
#else
#define MMPI_Comm int
#define MMPI_COMM_WORLD 0
#endif 

#include "Tessellation.hh"
#include <string>
#include <float.h>
#include <map>

namespace polytope
{

namespace Silo
{

// Helper function for finding available cycles.
std::vector<int> findAvailableCycles(const std::string& prefix,
                                     const std::string& directory,
                                     MMPI_Comm comm);

}

//! \class SiloReader
//! This class provides a static interface for reading Silo files 
//! containing tessellations made by polytope.
template <int Dimension, typename RealType>
class SiloReader
{
  // No general recipe
};

//! Partial specialization for 2D tessellations.
template <typename RealType>
class SiloReader<2, RealType>
{
  public:

  //! Returns a list of cycle numbers for Silo files dumped by a SiloWriter
  //! with the given prefix, in the given directory. If the directory is 
  //! omitted, its name is generated automatically from the prefix.
  static std::vector<int> availableCycles(const std::string& filePrefix,
                                          const std::string& directory = "",
                                          MMPI_Comm comm = MMPI_COMM_WORLD)
  {
    return Silo::findAvailableCycles(filePrefix, directory, comm);
  }

  //! Read an arbitrary polygonal mesh and an associated set of 
  //! fields and tags from a SILO file in the given directory.
  //! \param fields A map that will store arrays of field data read in from 
  //!               the file. If \a fields contains keys, only those fields
  //!               with those keys will be read from the file, and an error 
  //!               will occur if any of the keys are not found. If \a fields 
  //!               is empty, all data will be read in from the file.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void read(Tessellation<2, RealType>& mesh, 
                   std::map<std::string, std::vector<RealType> >& fields,
                   std::map<std::string, std::vector<int> >& nodeTags,
                   std::map<std::string, std::vector<int> >& edgeTags,
                   std::map<std::string, std::vector<int> >& faceTags,
                   std::map<std::string, std::vector<int> >& cellTags,
                   const std::string& filePrefix,
                   const std::string& directory,
                   int cycle,
                   RealType& time,
                   MMPI_Comm comm = MMPI_COMM_WORLD,
                   int numFiles = -1,
                   int mpiTag = 0);

  //! Read an arbitrary polygonal mesh and an associated set of 
  //! fields from a SILO file in the given directory.
  //! \param fields A map that will store arrays of field data read in from 
  //!               the file. If \a fields contains keys, only those fields
  //!               with those keys will be read from the file, and an error 
  //!               will occur if any of the keys are not found. If \a fields 
  //!               is empty, all data will be read in from the file.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void read(Tessellation<2, RealType>& mesh, 
                   std::map<std::string, std::vector<RealType> >& fields,
                   const std::string& filePrefix,
                   const std::string& directory,
                   int cycle,
                   RealType& time,
                   MMPI_Comm comm = MMPI_COMM_WORLD,
                   int numFiles = -1,
                   int mpiTag = 0)
  {
    std::map<std::string, std::vector<int> > nTags, eTags, fTags, cTags;
    read(mesh, fields, nTags, eTags, fTags, cTags, filePrefix, directory,
         cycle, time, comm, numFiles, mpiTag);
  }

  //! Read an arbitrary polygonal mesh and an associated set of 
  //! fields from a SILO file. This version generates a 
  //! directory name automatically. For parallel runs, the directory 
  //! name is filePrefix-nproc. For serial runs, the directory is 
  //! the current working directory.
  //! \param fields A map that will store arrays of field data read in from 
  //!               the file. If \a fields contains keys, only those fields
  //!               with those keys will be read from the file, and an error 
  //!               will occur if any of the keys are not found. If \a fields 
  //!               is empty, all data will be read in from the file.
  //!                 is set to -1, one file will be written for each process.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void read(Tessellation<2, RealType>& mesh, 
                   std::map<std::string, std::vector<RealType> >& fields,
                   const std::string& filePrefix,
                   int cycle,
                   RealType& time,
                   MMPI_Comm comm = MMPI_COMM_WORLD,
                   int numFiles = -1,
                   int mpiTag = 0)
  {
    read(mesh, fields, filePrefix, "", cycle, time, comm, numFiles, mpiTag);
  }

  //! This version of read omits the cycle and time arguments.
  static void read(Tessellation<2, RealType>& mesh, 
                   std::map<std::string, std::vector<RealType> >& fields,
                   const std::string& filePrefix,
                   const std::string& directory,
                   MMPI_Comm comm = MMPI_COMM_WORLD,
                   int numFiles = -1,
                   int mpiTag = 0)
  {
    RealType time;
    read(mesh, fields, filePrefix, directory, -1, time,
         comm, numFiles, mpiTag);
  }

  //! This version of read omits the cycle and time arguments and 
  //! automatically generates the directory name from the file prefix.
  static void read(Tessellation<2, RealType>& mesh, 
                   std::map<std::string, std::vector<RealType> >& fields,
                   const std::string& filePrefix,
                   MMPI_Comm comm = MMPI_COMM_WORLD,
                   int numFiles = -1,
                   int mpiTag = 0)
  {
    RealType time;
    read(mesh, fields, filePrefix, "", -1, time,
         comm, numFiles, mpiTag);
  }

};

//! Partial specialization for 3D tessellations.
template <typename RealType>
class SiloReader<3, RealType>
{
  public:

  //! Returns a list of cycle numbers for Silo files dumped by a SiloWriter
  //! with the given prefix, in the given directory. If the directory is 
  //! omitted, its name is generated automatically from the prefix.
  static std::vector<int> availableCycles(const std::string& filePrefix,
                                          const std::string& directory = "",
                                          MMPI_Comm comm = MMPI_COMM_WORLD)
  {
    return Silo::findAvailableCycles(filePrefix, directory, comm);
  }

  //! Read an arbitrary polyhedral mesh and an associated set of 
  //! fields and tags from a SILO file in the given directory.
  //! \param fields A map that will store arrays of field data read in from 
  //!               the file. If \a fields contains keys, only those fields
  //!               with those keys will be read from the file, and an error 
  //!               will occur if any of the keys are not found. If \a fields 
  //!               is empty, all data will be read in from the file.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void read(Tessellation<3, RealType>& mesh, 
                   std::map<std::string, std::vector<RealType> >& fields,
                   std::map<std::string, std::vector<int> >& nodeTags,
                   std::map<std::string, std::vector<int> >& edgeTags,
                   std::map<std::string, std::vector<int> >& faceTags,
                   std::map<std::string, std::vector<int> >& cellTags,
                   const std::string& filePrefix,
                   const std::string& directory,
                   int cycle,
                   RealType& time,
                   MMPI_Comm comm = MMPI_COMM_WORLD,
                   int numFiles = -1,
                   int mpiTag = 0);

  //! Read an arbitrary polyhedral mesh and an associated set of 
  //! fields from a SILO file in the given directory.
  //! \param fields A map that will store arrays of field data read in from 
  //!               the file. If \a fields contains keys, only those fields
  //!               with those keys will be read from the file, and an error 
  //!               will occur if any of the keys are not found. If \a fields 
  //!               is empty, all data will be read in from the file.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void read(Tessellation<3, RealType>& mesh, 
                   std::map<std::string, std::vector<RealType> >& fields,
                   const std::string& filePrefix,
                   const std::string& directory,
                   int cycle,
                   RealType& time,
                   MMPI_Comm comm = MMPI_COMM_WORLD,
                   int numFiles = -1,
                   int mpiTag = 0)
  {
    std::map<std::string, std::vector<int> > nTags, eTags, fTags, cTags;
    read(mesh, fields, nTags, eTags, fTags, cTags, filePrefix, directory,
         cycle, time, comm, numFiles, mpiTag);
  }

  //! Read an arbitrary polyhedral mesh and an associated set of 
  //! fields from a SILO file. This version generates a 
  //! directory name automatically. For parallel runs, the directory 
  //! name is filePrefix-nproc. For serial runs, the directory is 
  //! the current working directory.
  //! \param fields A map that will store arrays of field data read in from 
  //!               the file. If \a fields contains keys, only those fields
  //!               with those keys will be read from the file, and an error 
  //!               will occur if any of the keys are not found. If \a fields 
  //!               is empty, all data will be read in from the file.
  //! \param numFiles The number of files that will be written. If this 
  //!                 is set to -1, one file will be written for each process.
  static void read(Tessellation<3, RealType>& mesh, 
                   std::map<std::string, std::vector<RealType> >& fields,
                   const std::string& filePrefix,
                   int cycle,
                   RealType& time,
                   MMPI_Comm comm = MMPI_COMM_WORLD,
                   int numFiles = -1,
                   int mpiTag = 0)
  {
    read(mesh, fields, filePrefix, "", cycle, time, comm, numFiles, mpiTag);
  }

  //! This version of read omits the cycle and time arguments.
  static void read(Tessellation<3, RealType>& mesh, 
                   std::map<std::string, std::vector<RealType> >& fields,
                   const std::string& filePrefix,
                   const std::string& directory,
                   MMPI_Comm comm = MMPI_COMM_WORLD,
                   int numFiles = -1,
                   int mpiTag = 0)
  {
    RealType time;
    read(mesh, fields, filePrefix, directory, -1, time,
         comm, numFiles, mpiTag);
  }

  //! This version of read omits the cycle and time arguments and 
  //! automatically generates the directory name from the file prefix.
  static void read(Tessellation<3, RealType>& mesh, 
                   std::map<std::string, std::vector<RealType> >& fields,
                   const std::string& filePrefix,
                   MMPI_Comm comm = MMPI_COMM_WORLD,
                   int numFiles = -1,
                   int mpiTag = 0)
  {
    RealType time;
    read(mesh, fields, filePrefix, "", -1, time,
         comm, numFiles, mpiTag);
  }

};

} // end namespace

#endif
