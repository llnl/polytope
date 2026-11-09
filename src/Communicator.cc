#include "Communicator.hh"
#include "polytope_internal.hh"

#include <algorithm>
#include <cstdlib>
#include <iostream>

namespace {

#ifdef POLYTOPE_ENABLE_MPI
bool mpiIsActive() {
  int initialized = 0;
  if (MPI_Initialized(&initialized) != MPI_SUCCESS || !initialized) return false;

  int finalized = 0;
  return MPI_Finalized(&finalized) == MPI_SUCCESS && !finalized;
}

[[noreturn]] void handleMpiError(const char* operation, const int error) {
  const bool mpiActive = mpiIsActive();
  char errorString[MPI_MAX_ERROR_STRING + 1] = {};
  int errorLength = 0;
  if (!mpiActive || MPI_Error_string(error, errorString, &errorLength) != MPI_SUCCESS) {
    std::cerr << operation << " failed with MPI error code " << error << std::endl;
  } else {
    errorString[std::min(errorLength, MPI_MAX_ERROR_STRING)] = '\0';
    std::cerr << operation << " failed: " << errorString << std::endl;
  }
  std::cerr.flush();

  if (mpiActive) MPI_Abort(polytope::Communicator::communicator(), error);
  std::abort();
}

void checkMpi(const int error, const char* operation) {
  if (error != MPI_SUCCESS) handleMpiError(operation, error);
}
#endif

}

namespace polytope {

//------------------------------------------------------------------------------
// Default constructor (private).
//------------------------------------------------------------------------------
Communicator::Communicator() {
#ifdef POLYTOPE_ENABLE_MPI
  // Use world by default but can be changed
  mCommunicator = MPI_COMM_WORLD;
#else
  mCommunicator = 0;
#endif
}

//------------------------------------------------------------------------------
// Destructor (private).
//------------------------------------------------------------------------------
Communicator::~Communicator() {
}

//------------------------------------------------------------------------------
// Instance
//------------------------------------------------------------------------------
Communicator& Communicator::instance() {
  static Communicator theInstance;
  return theInstance;
}

//------------------------------------------------------------------------------
// Public routines
//------------------------------------------------------------------------------

void Communicator::init() {
#ifdef POLYTOPE_ENABLE_MPI
  int isInit;
  checkMpi(MPI_Initialized(&isInit), "MPI_Initialized");
  if (!isInit) {
    checkMpi(MPI_Init(nullptr, nullptr), "MPI_Init");
  }
#endif
}

void Communicator::init(int argc, char** argv) {
#ifdef POLYTOPE_ENABLE_MPI
  int isInit;
  checkMpi(MPI_Initialized(&isInit), "MPI_Initialized");
  if (!isInit) {
    checkMpi(MPI_Init(&argc, &argv), "MPI_Init");
  }
#else
  POLY_CONTRACT_VAR(argc);
  POLY_CONTRACT_VAR(argv);
#endif
}

MPI_Comm* Communicator::comm_ptr() {
#ifdef POLYTOPE_ENABLE_MPI
  return &(instance().mCommunicator);
#else
  return nullptr;
#endif
}

void Communicator::finalize() {
#ifdef POLYTOPE_ENABLE_MPI
  if (mpiIsActive()) checkMpi(MPI_Finalize(), "MPI_Finalize");
#endif
}

int Communicator::getRank() {
  int sRank = 0;
#ifdef POLYTOPE_ENABLE_MPI
  if (mpiIsActive()) {
    checkMpi(MPI_Comm_rank(communicator(), &sRank), "MPI_Comm_rank");
  }
#endif
  return sRank;
}

int Communicator::getNRanks() {
  int nRanks = 1;
#ifdef POLYTOPE_ENABLE_MPI
  if (mpiIsActive()) {
    checkMpi(MPI_Comm_size(communicator(), &nRanks), "MPI_Comm_size");
  }
#endif
  return nRanks;
}

void Communicator::Barrier() {
#ifdef POLYTOPE_ENABLE_MPI
  if (mpiIsActive()) {
    checkMpi(MPI_Barrier(communicator()), "MPI_Barrier");
  }
#endif
}

[[noreturn]] void Communicator::abort() {
#ifdef POLYTOPE_ENABLE_MPI
  if (mpiIsActive()) {
    MPI_Abort(communicator(), 1);
  }
#endif
  std::abort();
}

int Communicator::getRoot() {
  return instance().m_root;
}

void Communicator::setRoot(const int root) {
  instance().m_root = root;
}

}
