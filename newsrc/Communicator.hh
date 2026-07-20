// Communicator
// A singleton object which holds the MPI communicator. Taken from Spheral.
//
#ifndef __Polytope_Communicator__
#define __Polytope_Communicator__

#include "polytope.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

namespace polytope {

class Communicator {

public:
  //------------------------===== Public Interface =====-----------------------//
  // Get the instance.
  static Communicator& instance();
  // Run MPI_Init
  static void init();
  static void init(int argc, char** argv);
  // Run MPI_Finalize
  static void finalize();

  // Access the communicator.
  static MPI_Comm& communicator() { return instance().mCommunicator; }
  static void communicator(MPI_Comm& comm) { instance().mCommunicator = comm; }
  static MPI_Comm* comm_ptr();
  static int getRank();
  static int getNProcs();
  static void Barrier();
  static void haltAll();
  static int getRoot();
  static void setRoot(const int root);

private:
  //------------------------===== Private Interface =====----------------------//
  MPI_Comm mCommunicator;
  int m_root = 0;

  // No public constructors, destructor, or assignment.
  Communicator();
  Communicator(const Communicator&);
  Communicator& operator=(const Communicator&);
  ~Communicator();
};

}

#endif
