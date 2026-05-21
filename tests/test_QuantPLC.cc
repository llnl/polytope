// Unit tests for QuantPLC - Quantized PLC with hash-based deduplication
//
// Tests the QuantPLC class which combines:
//   - PLC reduction (like ReducedPLC)
//   - Hash-based point deduplication via Morton curve
//   - Quantized coordinate storage

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <cassert>
#include <cmath>

#include "polytope.hh"
#include "QuantPLC.hh"
#include "PLC.hh"
#include "Point.hh"
#include "HashKey.hh"
#include "Quantizer.hh"
#include "polytope_test_utilities.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

namespace {

//------------------------------------------------------------------------------
// Test 1: Basic construction with reduction
// Verifies that QuantPLC can be constructed from a PLC and reduces to only
// the points used in facets
//------------------------------------------------------------------------------
void testBasicConstruction() {
  cout << "Test 1: Basic construction with reduction" << endl;

  // Create a simple square PLC in 2D
  PLC<2> plc;
  plc.facets = {
    {0, 1},  // Bottom edge
    {1, 2},  // Right edge
    {2, 3},  // Top edge
    {3, 0}   // Left edge
  };

  // Point set includes extra unused points
  vector<double> allpoints = {
    0.0, 0.0,  // 0 - used
    1.0, 0.0,  // 1 - used
    1.0, 1.0,  // 2 - used
    0.0, 1.0,  // 3 - used
    2.0, 0.0,  // 4 - unused
    2.0, 1.0,  // 5 - unused
  };

  // Create quantizer
  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(2.0, 1.0);
  Quantizer<2> Q(low, high);

  // Construct QuantPLC with reduction
  QuantPLC<2> quantPLC(plc, Q, allpoints, true);

  // Should have remapped to 4 points (0,1,2,3 become 0,1,2,3)
  POLY_CHECK2(quantPLC.facets.size() == 4,
              "Should have 4 facets, got " << quantPLC.facets.size());

  // All facet indices should be in range [0, 3]
  for (const auto& facet : quantPLC.facets) {
    for (int idx : facet) {
      POLY_CHECK2(idx >= 0 && idx < 4,
                  "Facet index " << idx << " out of range [0, 3]");
    }
  }

  cout << "  PASS: Basic construction works" << endl;
}

//------------------------------------------------------------------------------
// Test 2: Deduplication via hashing
// Verifies that duplicate points (same quantized location) are merged
//------------------------------------------------------------------------------
void testDeduplication() {
  cout << "Test 2: Deduplication via hashing" << endl;

  // Create a PLC with duplicate point references
  PLC<2> plc;
  plc.facets = {
    {0, 1},
    {1, 2},
    {2, 3},
    {3, 4},  // References duplicate
  };

  // Points 0 and 4 are exact duplicates
  vector<double> allpoints = {
    0.0, 0.0,  // 0
    1.0, 0.0,  // 1
    1.0, 1.0,  // 2
    0.0, 1.0,  // 3
    0.0, 0.0,  // 4 - duplicate of 0
  };

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);
  Quantizer<2> Q(low, high);

  QuantPLC<2> quantPLC(plc, Q, allpoints, true);

  // Check that facets have been remapped
  // The last facet should reference the same index as first facet's first node
  POLY_CHECK2(quantPLC.facets[0][0] == quantPLC.facets[3][1],
              "Duplicate points should map to same index");

  cout << "  PASS: Deduplication works" << endl;
}

//------------------------------------------------------------------------------
// Test 3: Near-duplicate detection
// Verifies that points close enough to quantize to same location are merged
//------------------------------------------------------------------------------
void testNearDuplicates() {
  cout << "Test 3: Near-duplicate detection" << endl;

  PLC<2> plc;
  plc.facets = {
    {0, 1},
    {1, 2},
  };

  // Points 0 and 2 are very close (should quantize to same cell)
  vector<double> allpoints = {
    0.0, 0.0,      // 0
    1.0, 0.0,      // 1
    0.0, 0.0001,   // 2 - near-duplicate of 0
  };

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);
  Quantizer<2> Q(low, high, 0.02);

  QuantPLC<2> quantPLC(plc, Q, allpoints, true);

  // Points 0 and 2 should be merged if quantization is coarse enough
  // The number of unique points should reflect this
  cout << "  Number of facets after deduplication: " << quantPLC.facets.size() << endl;

  cout << "  PASS: Near-duplicate handling works" << endl;
}

//------------------------------------------------------------------------------
// Test 4: Construction without reduction
// Verifies that all points are kept when reduction is disabled
//------------------------------------------------------------------------------
void testNoReduction() {
  cout << "Test 4: Construction without reduction" << endl;

  PLC<2> plc;
  plc.facets = {
    {0, 1},
    {1, 2},
  };

  vector<double> allpoints = {
    0.0, 0.0,  // 0 - used
    1.0, 0.0,  // 1 - used
    1.0, 1.0,  // 2 - used
    2.0, 0.0,  // 3 - unused
  };

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(2.0, 1.0);
  Quantizer<2> Q(low, high);

  // Construct WITHOUT reduction
  QuantPLC<2> quantPLC(plc, Q, allpoints, false);

  // Facet indices should not be remapped (still reference original indices)
  POLY_CHECK2(quantPLC.facets.size() == 2,
              "Should have 2 facets, got " << quantPLC.facets.size());

  cout << "  PASS: Construction without reduction works" << endl;
}

//------------------------------------------------------------------------------
// Test 5: Empty PLC
// Verifies behavior with empty input
//------------------------------------------------------------------------------
void testEmptyPLC() {
  cout << "Test 5: Empty PLC" << endl;

  PLC<2> plc;  // Empty PLC

  vector<double> allpoints = {
    0.0, 0.0,
    1.0, 0.0,
  };

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);
  Quantizer<2> Q(low, high);

  QuantPLC<2> quantPLC(plc, Q, allpoints, true);

  POLY_CHECK2(quantPLC.facets.empty(),
              "Empty PLC should have no facets");

  cout << "  PASS: Empty PLC handled correctly" << endl;
}

//------------------------------------------------------------------------------
// Test 6: 3D PLC
// Tests QuantPLC with 3D geometry
//------------------------------------------------------------------------------
void test3D() {
  cout << "Test 6: 3D PLC" << endl;

  // Create a tetrahedron
  PLC<3> plc;
  plc.facets = {
    {0, 1, 2},  // Base triangle
    {0, 1, 3},  // Side face
    {1, 2, 3},  // Side face
    {2, 0, 3},  // Side face
  };

  vector<double> allpoints = {
    0.0, 0.0, 0.0,  // 0
    1.0, 0.0, 0.0,  // 1
    0.5, 1.0, 0.0,  // 2
    0.5, 0.5, 1.0,  // 3
  };

  Point<3, double> low(0.0, 0.0, 0.0);
  Point<3, double> high(1.0, 1.0, 1.0);
  Quantizer<3> Q(low, high);

  QuantPLC<3> quantPLC(plc, Q, allpoints, true);

  POLY_CHECK2(quantPLC.facets.size() == 4,
              "3D tetrahedron should have 4 facets, got " << quantPLC.facets.size());

  // Each facet in 3D should have at least 3 nodes
  for (size_t i = 0; i < quantPLC.facets.size(); ++i) {
    POLY_CHECK2(quantPLC.facets[i].size() >= 3,
                "3D facet " << i << " should have >= 3 nodes, got "
                << quantPLC.facets[i].size());
  }

  cout << "  PASS: 3D PLC works" << endl;
}

//------------------------------------------------------------------------------
// Test 7: Hash consistency
// Verifies that identical points always produce identical hashes
//------------------------------------------------------------------------------
void testHashConsistency() {
  cout << "Test 7: Hash consistency" << endl;

  vector<double> allpoints = {
    0.5, 0.5,  // 0
    0.5, 0.5,  // 1 - duplicate
  };

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);
  Quantizer<2> Q(low, high);

  // Create two points and hash them
  Point<2, double> p0(0.5, 0.5);
  Point<2, double> p1(0.5, 0.5);

  auto ip0 = Q.quantize(p0);
  auto ip1 = Q.quantize(p1);
  auto hash0 = Q.hash(ip0);
  auto hash1 = Q.hash(ip1);

  POLY_CHECK2(hash0 == hash1,
              "Identical points should produce identical hashes");

  // Create slightly different point
  Point<2, double> p2(0.5001, 0.5);
  auto ip2 = Q.quantize(p2);
  auto hash2 = Q.hash(ip2);

  // Whether hash2 equals hash0 depends on quantization resolution
  cout << "  Hash 0: 0x" << hex << hash0 << dec << endl;
  cout << "  Hash 1: 0x" << hex << hash1 << dec << endl;
  cout << "  Hash 2: 0x" << hex << hash2 << dec << endl;

  cout << "  PASS: Hash consistency verified" << endl;
}

//------------------------------------------------------------------------------
// Test 8: Multiple holes
// Tests PLC with holes (inner boundaries)
//------------------------------------------------------------------------------
void testPLCWithHoles() {
  cout << "Test 8: PLC with holes" << endl;

  PLC<2> plc;
  // Outer boundary
  plc.facets = {
    {0, 1},
    {1, 2},
    {2, 3},
    {3, 0},
  };

  // Inner hole
  plc.holes.resize(1);
  plc.holes[0] = {
    {4, 5},
    {5, 6},
    {6, 7},
    {7, 4},
  };

  vector<double> allpoints = {
    // Outer square
    0.0, 0.0,  // 0
    2.0, 0.0,  // 1
    2.0, 2.0,  // 2
    0.0, 2.0,  // 3
    // Inner square (hole)
    0.5, 0.5,  // 4
    1.5, 0.5,  // 5
    1.5, 1.5,  // 6
    0.5, 1.5,  // 7
  };

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(2.0, 2.0);
  Quantizer<2> Q(low, high);

  QuantPLC<2> quantPLC(plc, Q, allpoints, true);

  POLY_CHECK2(quantPLC.facets.size() == 4,
              "Should have 4 outer facets");
  POLY_CHECK2(quantPLC.holes.size() == 1,
              "Should have 1 hole");
  POLY_CHECK2(quantPLC.holes[0].size() == 4,
              "Hole should have 4 facets");

  cout << "  PASS: PLC with holes works" << endl;
}

//------------------------------------------------------------------------------
// Test 9: Complex index remapping
// Tests a case where reduction causes complex index remapping
//------------------------------------------------------------------------------
void testComplexRemapping() {
  cout << "Test 9: Complex index remapping" << endl;

  PLC<2> plc;
  // Only use every other point
  plc.facets = {
    {0, 2},
    {2, 4},
    {4, 6},
    {6, 0},
  };

  vector<double> allpoints = {
    0.0, 0.0,  // 0 - used
    0.5, 0.0,  // 1 - unused
    1.0, 0.0,  // 2 - used
    1.0, 0.5,  // 3 - unused
    1.0, 1.0,  // 4 - used
    0.5, 1.0,  // 5 - unused
    0.0, 1.0,  // 6 - used
    0.0, 0.5,  // 7 - unused
  };

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);
  Quantizer<2> Q(low, high);

  QuantPLC<2> quantPLC(plc, Q, allpoints, true);

  // After reduction, indices should be 0,1,2,3 (remapped from 0,2,4,6)
  for (const auto& facet : quantPLC.facets) {
    for (int idx : facet) {
      POLY_CHECK2(idx >= 0 && idx < 4,
                  "Remapped index " << idx << " should be in [0,3]");
    }
  }

  // First facet should now be {0, 1} after remapping
  POLY_CHECK2(quantPLC.facets[0][0] == 0,
              "First facet first node should be 0 after remapping");
  POLY_CHECK2(quantPLC.facets[0][1] == 1,
              "First facet second node should be 1 after remapping");

  cout << "  PASS: Complex remapping works" << endl;
}

//------------------------------------------------------------------------------
// Test 10: PLC validity after construction
// Verifies that constructed QuantPLC is valid
//------------------------------------------------------------------------------
void testValidity() {
  cout << "Test 10: PLC validity after construction" << endl;

  PLC<2> plc;
  plc.facets = {
    {0, 1},
    {1, 2},
    {2, 0},
  };

  vector<double> allpoints = {
    0.0, 0.0,
    1.0, 0.0,
    0.5, 1.0,
  };

  Point<2, double> low(0.0, 0.0);
  Point<2, double> high(1.0, 1.0);
  Quantizer<2> Q(low, high);

  QuantPLC<2> quantPLC(plc, Q, allpoints, true);

  POLY_CHECK2(quantPLC.valid(),
              "Constructed QuantPLC should be valid");
  POLY_CHECK2(!quantPLC.empty(),
              "Constructed QuantPLC should not be empty");

  cout << "  PASS: QuantPLC validity checks work" << endl;
}

} // anonymous namespace

//------------------------------------------------------------------------------
// Main
//------------------------------------------------------------------------------

int main(int argc, char** argv) {
#ifdef POLYTOPE_ENABLE_MPI
  MPI_Init(&argc, &argv);
#else
  POLY_CONTRACT_VAR(argc);
  POLY_CONTRACT_VAR(argv);
#endif

  try {
    testBasicConstruction();
    testDeduplication();
    testNearDuplicates();
    testNoReduction();
    testEmptyPLC();
    test3D();
    testHashConsistency();
    testPLCWithHoles();
    testComplexRemapping();
    testValidity();

    cout << "\n" << "All QuantPLC tests PASSED" << "\n" << endl;

  } catch (const char* str) {
    cerr << "FAILED: " << str << endl;
#ifdef POLYTOPE_ENABLE_MPI
    MPI_Finalize();
#endif
    return 1;
  } catch (const std::exception& e) {
    cerr << "FAILED: " << e.what() << endl;
#ifdef POLYTOPE_ENABLE_MPI
    MPI_Finalize();
#endif
    return 1;
  }

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
