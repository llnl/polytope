#include <iostream>
#include <cstdlib>

#include "polytope_internal.hh"

#ifdef POLYTOPE_ENABLE_MPI
// extern "C" {
#include <mpi.h>
// }
#endif

namespace polytope {

void internal_abort() {
  std::cout.flush();
  std::cerr.flush();

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Abort(MPI_COMM_WORLD, -1);
#else
  abort();
#endif
}

}
