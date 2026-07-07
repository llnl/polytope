#include <iostream>
#include <cstdlib>

#include "polytope_internal.hh"
#include "Communicator.hh"

namespace polytope {

void internal_abort() {

#ifdef POLYTOPE_ENABLE_MPI
  const int rank = Communicator::getRank();
  if (rank == 0) {
    std::cout.flush();
    std::cerr.flush();
  }
  Communicator::haltAll();
  abort();
#else
  std::cout.flush();
  std::cerr.flush();
  abort();
#endif
}

}
