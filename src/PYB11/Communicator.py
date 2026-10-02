from PYB11Generator import *

@PYB11singleton
class Communicator:
    """
Singleton for doing MPI communication.
The init and finalize are called in the polytope module itself and are not provided here.
"""

    @PYB11static
    @PYB11returnpolicy("reference")
    def instance(self):
        return "Communicator&"

    @PYB11static
    def getRank(self):
        return "int"

    @PYB11static
    def getNRanks(self):
        return "int"

    @PYB11static
    def getRoot(self):
        return "int"

    @PYB11static
    def setRoot(self,
                root="const int"):
        return "void"

    @PYB11static
    def Barrier(self):
        return "void"

    @PYB11static
    def abort(self):
        return "void"
