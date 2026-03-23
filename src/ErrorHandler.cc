#include "polytope.hh"

#include <iostream>
#include <cstdlib>

#ifdef POLYTOPE_ENABLE_MPI
// extern "C" {
#include "mpi.h"
// }
#endif

namespace polytope
{

namespace {

// Error handler.
static void (*errorHandler)(const std::string&, int) = 0;

//------------------------------------------------------------------------
// Default error handler.
//------------------------------------------------------------------------
void defaultErrorHandler(const std::string& message, int status)
{
#ifdef POLYTOPE_ENABLE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
  {
    std::cout << "Error: " << message << std::endl;
    std::cout << "encountered in polytope library. Exiting with status " << status << std::endl;
  }
  MPI_Abort(MPI_COMM_WORLD, status);
#else
  std::cout << "Error: " << message << std::endl;
  std::cout << "encountered in polytope library. Exiting with status " << status << std::endl;
  exit(status);
#endif
}
//------------------------------------------------------------------------

} // end anonymous namespace

//------------------------------------------------------------------------
void 
error(const std::string& message, int status)
{
  // Set the default error handler if no handler is set.
  if (errorHandler == 0)
    errorHandler = defaultErrorHandler;

  // Call it.
  errorHandler(message, status);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void 
setErrorHandler(void (*handler)(const std::string&, int))
{
  errorHandler = handler;
}
//------------------------------------------------------------------------

}
