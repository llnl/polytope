#include <iostream>
#include <cstdlib>

#include "polytope_internal.hh"
#include "Communicator.hh"

namespace polytope {

[[noreturn]] void internal_abort() {

#ifdef POLYTOPE_ENABLE_MPI
  std::cout.flush();
  std::cerr.flush();
  Communicator::abort();
#else
  std::cout.flush();
  std::cerr.flush();
  std::abort();
#endif
}

}
