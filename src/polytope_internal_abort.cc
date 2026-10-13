#include <iostream>
#include <cstdlib>

#include "polytope_internal.hh"
#include "Communicator.hh"

namespace polytope {

void internal_abort() {

#ifdef POLYTOPE_ENABLE_MPI
  const int rank = Communicator::getRank();
  if (rank == Communicator::getRoot()) {
    std::cout.flush();
    std::cerr.flush();
  }
  Communicator::abort();
#else
  std::cout.flush();
  std::cerr.flush();
  abort();
#endif
}

}
