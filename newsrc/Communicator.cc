#include "Communicator.hh"
#include "polytope_internal.hh"

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
  MPI_Initialized(&isInit);
  if (!isInit) {
    MPI_Init(nullptr, nullptr);
  }
#endif
}

void Communicator::init(int argc, char** argv) {
#ifdef POLYTOPE_ENABLE_MPI
  int isInit;
  MPI_Initialized(&isInit);
  if (!isInit) {
    MPI_Init(&argc, &argv);
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
  int finalized = 0;
  MPI_Finalized(&finalized);
  if (finalized == 0) {
    int finalize = MPI_Finalize();
    if (finalize != 0) {
      char string[MPI_MAX_ERROR_STRING];
      int resultlen = 0;
      MPI_Error_string(finalize, string, &resultlen);
      POLY_VERIFY2(finalize == 0, string);
    }
  }
#endif
}

int Communicator::getRank() {
  int sRank = 0;
#ifdef POLYTOPE_ENABLE_MPI
  int isInit;
  MPI_Initialized(&isInit);
  if (isInit) {
    MPI_Comm_rank(communicator(), &sRank);
  }
#endif
  return sRank;
}

int Communicator::getNProcs() {
  int nRanks = 1;
#ifdef POLYTOPE_ENABLE_MPI
  int isInit;
  MPI_Initialized(&isInit);
  if (isInit) {
    MPI_Comm_size(communicator(), &nRanks);
  }
#endif
  return nRanks;
}

void Communicator::Barrier() {
#ifdef POLYTOPE_ENABLE_MPI
  int isInit;
  MPI_Initialized(&isInit);
  if (isInit) {
    MPI_Barrier(communicator());
  }
#endif
}

void Communicator::abort() {
#ifdef POLYTOPE_ENABLE_MPI
  int isInit;
  MPI_Initialized(&isInit);
  if (isInit) {
    MPI_Abort(communicator(), 1);
  }
#endif
}

int Communicator::getRoot() {
  return instance().m_root;
}

void Communicator::setRoot(const int root) {
  instance().m_root = root;
}

}
