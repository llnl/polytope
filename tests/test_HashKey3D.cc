// Unit tests for HashKey3D Morton encoding/decoding (128-bit)
//
// NOTE: HashKey3D only supports non-negative coordinates.
// The implementation uses unsigned integer bit manipulation internally.

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <cassert>
#include <cstdint>
#include <algorithm>
#include <iomanip>

#include "polytope.hh"
#include "HashKey.hh"
#include "Point.hh"
#include "polytope_test_utilities.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

namespace {

//------------------------------------------------------------------------------
// Type alias for HashKey3D
//------------------------------------------------------------------------------
using HashKey3D = HashKey<3, double>;

//------------------------------------------------------------------------------
// Test that hash/unhash are inverses
//------------------------------------------------------------------------------
void testRoundTrip(const Point3<uint64_t>& p, const string& label) {
  auto key = HashKey3D::hash(p);
  Point3<uint64_t> p2 = HashKey3D::unhash(key);

  POLY_CHECK2(p.x == p2.x && p.y == p2.y && p.z == p2.z,
              label << ": Round-trip failed for (" << p.x << ", " << p.y << ", " << p.z
              << ") -> (" << p2.x << ", " << p2.y << ", " << p2.z << ")");
}

//------------------------------------------------------------------------------
// Test that identical points produce identical hashes
//------------------------------------------------------------------------------
void testIdenticalPoints(const Point3<uint64_t>& p) {
  auto key1 = HashKey3D::hash(p);
  auto key2 = HashKey3D::hash(p);

  POLY_CHECK2(key1 == key2,
              "Identical points should hash identically: (" << p.x << ", " << p.y
              << ", " << p.z << ")");
}

//------------------------------------------------------------------------------
// Test that different points produce different hashes
//------------------------------------------------------------------------------
void testDistinctPoints(const Point3<uint64_t>& p1, const Point3<uint64_t>& p2) {
  if (p1.x == p2.x && p1.y == p2.y && p1.z == p2.z) return;

  auto key1 = HashKey3D::hash(p1);
  auto key2 = HashKey3D::hash(p2);

  POLY_CHECK2(key1 != key2,
              "Distinct points should hash differently: (" << p1.x << ", " << p1.y
              << ", " << p1.z << ") vs (" << p2.x << ", " << p2.y << ", " << p2.z << ")");
}

//------------------------------------------------------------------------------
// Test uniqueness - N distinct points produce N distinct hashes
//------------------------------------------------------------------------------
void testUniqueness(const vector<Point3<uint64_t>>& points) {
  map<unsigned __int128, vector<unsigned>> collisions;

  for (unsigned i = 0; i < points.size(); ++i) {
    auto key = HashKey3D::hash(points[i]);
    collisions[key].push_back(i);
  }

  // Report any collisions
  unsigned collisionCount = 0;
  for (const auto& kv : collisions) {
    if (kv.second.size() > 1) {
      collisionCount++;
      uint64_t lo = (uint64_t)kv.first;
      uint64_t hi = (uint64_t)(kv.first >> 64);
      cout << "  Collision at hash (" << hex << lo << ", " << hi << dec << "): ";
      for (unsigned idx : kv.second) {
        cout << "(" << points[idx].x << "," << points[idx].y << ","
             << points[idx].z << ") ";
      }
      cout << endl;
    }
  }

  POLY_CHECK2(collisionCount == 0,
              "Found " << collisionCount << " hash collisions among "
              << points.size() << " distinct points");
}

//------------------------------------------------------------------------------
// Test spatial locality for 3D Morton curves
// Uses a grid of points to ensure meaningful locality testing
//------------------------------------------------------------------------------
void testSpatialLocality() {
  cout << "Testing spatial locality with 3D grid..." << endl;

  // Create a regular 3D grid of points (16x16x16 = 4096 points)
  const unsigned gridSize = 16;
  vector<Point3<uint64_t>> points;
  points.reserve(gridSize * gridSize * gridSize);

  for (unsigned i = 0; i < gridSize; ++i) {
    for (unsigned j = 0; j < gridSize; ++j) {
      for (unsigned k = 0; k < gridSize; ++k) {
        points.push_back(Point3<uint64_t>(uint64_t(i), uint64_t(j), uint64_t(k)));
      }
    }
  }

  // Hash and sort by Morton code
  vector<pair<unsigned __int128, unsigned>> hashed;
  for (unsigned i = 0; i < points.size(); ++i) {
    auto key = HashKey3D::hash(points[i]);
    hashed.push_back(make_pair(key, i));
  }

  // Sort by hash value
  sort(hashed.begin(), hashed.end());

  // Compute average distance and locality
  double totalDist = 0.0;
  unsigned nearbyCount = 0;
  for (unsigned i = 0; i + 1 < hashed.size(); ++i) {
    unsigned idx1 = hashed[i].second;
    unsigned idx2 = hashed[i+1].second;

    double dx = double(points[idx1].x > points[idx2].x ?
                       points[idx1].x - points[idx2].x : points[idx2].x - points[idx1].x);
    double dy = double(points[idx1].y > points[idx2].y ?
                       points[idx1].y - points[idx2].y : points[idx2].y - points[idx1].y);
    double dz = double(points[idx1].z > points[idx2].z ?
                       points[idx1].z - points[idx2].z : points[idx2].z - points[idx1].z);
    double dist = dx + dy + dz; // Manhattan distance
    totalDist += dist;

    if (dist <= 6.0) nearbyCount++;
  }

  double avgDist = totalDist / (hashed.size() - 1);
  double localityRatio = double(nearbyCount) / (hashed.size() - 1);

  cout << "  Average Manhattan distance between consecutive hashes: " << avgDist << endl;
  cout << "  Spatial locality ratio (dist <= 6): " << localityRatio << " ("
       << nearbyCount << "/" << (hashed.size() - 1) << ")" << endl;

  // Morton curves should have good locality
  // For a 16x16x16 grid, random pairing would give avg distance ~24
  POLY_CHECK2(avgDist < 15.0,
              "Average distance too high: " << avgDist << " (Morton locality poor)");
  POLY_CHECK2(localityRatio > 0.35,
              "Spatial locality ratio too low: " << localityRatio);
}

//------------------------------------------------------------------------------
// Test edge cases for 3D
//------------------------------------------------------------------------------
void testEdgeCases() {
  cout << "Testing edge cases..." << endl;

  // Origin
  testRoundTrip(Point3<uint64_t>(0, 0, 0), "Origin");

  // Max values for coordinate type (42 bits per dimension)
  const uint64_t maxVal = HashKey3D::coordMax();
  testRoundTrip(Point3<uint64_t>(maxVal, maxVal, maxVal), "Max coordinates");
  testRoundTrip(Point3<uint64_t>(maxVal, 0, 0), "Max X");
  testRoundTrip(Point3<uint64_t>(0, maxVal, 0), "Max Y");
  testRoundTrip(Point3<uint64_t>(0, 0, maxVal), "Max Z");
  testRoundTrip(Point3<uint64_t>(maxVal, maxVal, 0), "Max X,Y");
  testRoundTrip(Point3<uint64_t>(maxVal, 0, maxVal), "Max X,Z");
  testRoundTrip(Point3<uint64_t>(0, maxVal, maxVal), "Max Y,Z");

  // Powers of 2
  for (unsigned bit = 0; bit < min(20u, HashKey3D::num1DBits()); ++bit) {
    uint64_t val = uint64_t(1) << bit;
    testRoundTrip(Point3<uint64_t>(val, 0, 0), "Power of 2 in X");
    testRoundTrip(Point3<uint64_t>(0, val, 0), "Power of 2 in Y");
    testRoundTrip(Point3<uint64_t>(0, 0, val), "Power of 2 in Z");
    testRoundTrip(Point3<uint64_t>(val, val, val), "Power of 2 in all");
  }

  // Small integer grid
  for (uint64_t x = 0; x < 8; ++x) {
    for (uint64_t y = 0; y < 8; ++y) {
      for (uint64_t z = 0; z < 8; ++z) {
        testRoundTrip(Point3<uint64_t>(x, y, z), "Small integers");
        testIdenticalPoints(Point3<uint64_t>(x, y, z));

        if (x + 1 < 8) {
          testDistinctPoints(Point3<uint64_t>(x, y, z),
                           Point3<uint64_t>(x+1, y, z));
        }
        if (y + 1 < 8) {
          testDistinctPoints(Point3<uint64_t>(x, y, z),
                           Point3<uint64_t>(x, y+1, z));
        }
        if (z + 1 < 8) {
          testDistinctPoints(Point3<uint64_t>(x, y, z),
                           Point3<uint64_t>(x, y, z+1));
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// Test the 128-bit representation structure
//------------------------------------------------------------------------------
void test128BitStructure() {
  cout << "\n=== Testing 128-bit structure ===" << endl;

  // Test that we can represent 42 bits per dimension = 126 bits total
  Point3<uint64_t> p(0, 0, 0);

  // Set high bit in each dimension
  p.x = uint64_t(1) << 41;
  p.y = uint64_t(1) << 41;
  p.z = uint64_t(1) << 41;

  auto key = HashKey3D::hash(p);
  Point3<uint64_t> p2 = HashKey3D::unhash(key);

  POLY_CHECK2(p.x == p2.x && p.y == p2.y && p.z == p2.z,
              "Failed to represent 42-bit values: (" << p.x << ", " << p.y
              << ", " << p.z << ") -> (" << p2.x << ", " << p2.y << ", " << p2.z << ")");

  // Test that hi bits are used
  uint64_t lo = (uint64_t)key;
  uint64_t hi = (uint64_t)(key >> 64);
  cout << "  Hash of (2^41, 2^41, 2^41): lo=" << hex << lo
       << " hi=" << hi << dec << endl;
  POLY_CHECK2(hi != 0, "High 64 bits should be used for 42-bit coordinates");

  cout << "128-bit structure tests passed!" << endl;
}

} // anonymous namespace

//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Init(&argc, &argv);
#endif

  //----------------------------------------------------------------------------
  // Test 128-bit structure
  //----------------------------------------------------------------------------
  test128BitStructure();

  //----------------------------------------------------------------------------
  // Test 64-bit unsigned integers (using 42 bits per dimension)
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Testing HashKey3D with uint64_t ===" << endl;

    // Edge cases
    testEdgeCases();

    // Random points in 42-bit range
    cout << "Testing random points (42-bit range)..." << endl;
    const unsigned n = 800;
    vector<Point3<uint64_t>> points;
    for (unsigned i = 0; i < n; ++i) {
      uint64_t x = uint64_t(random01() * (1ULL << 40));
      uint64_t y = uint64_t(random01() * (1ULL << 40));
      uint64_t z = uint64_t(random01() * (1ULL << 40));
      points.push_back(Point3<uint64_t>(x, y, z));
      testRoundTrip(points.back(), "Random point");
    }

    testUniqueness(points);
    testSpatialLocality();
  }

  //----------------------------------------------------------------------------
  // Test with smaller range (20-bit values)
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Testing HashKey3D with 20-bit range ===" << endl;

    const unsigned n = 800;
    vector<Point3<uint64_t>> points;
    for (unsigned i = 0; i < n; ++i) {
      uint64_t x = uint64_t(random01() * (1ULL << 20));
      uint64_t y = uint64_t(random01() * (1ULL << 20));
      uint64_t z = uint64_t(random01() * (1ULL << 20));
      points.push_back(Point3<uint64_t>(x, y, z));
      testRoundTrip(points.back(), "Random 20-bit point");
    }

    testUniqueness(points);
    testSpatialLocality();
  }

  //----------------------------------------------------------------------------
  // Stress test
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Stress test ===" << endl;
    const unsigned nOps = 5000;
    cout << "Performing " << nOps << " hash/unhash operations..." << endl;

    for (unsigned i = 0; i < nOps; ++i) {
      uint64_t x = uint64_t(random01() * (1ULL << 40));
      uint64_t y = uint64_t(random01() * (1ULL << 40));
      uint64_t z = uint64_t(random01() * (1ULL << 40));

      auto key = HashKey3D::hash(Point3<uint64_t>(x, y, z));
      Point3<uint64_t> recovered = HashKey3D::unhash(key);

      POLY_CHECK(recovered.x == x && recovered.y == y && recovered.z == z);
    }
    cout << "Stress test passed!" << endl;
  }

  //----------------------------------------------------------------------------
  // Operator tests
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Testing operators ===" << endl;

    auto key1 = HashKey3D::hash(Point3<uint64_t>(42, 73, 101));
    auto key2 = HashKey3D::hash(Point3<uint64_t>(42, 73, 101));
    auto key3 = HashKey3D::hash(Point3<uint64_t>(101, 73, 42));

    // Equality
    POLY_CHECK(key1 == key2);
    POLY_CHECK(!(key1 != key2));

    // Inequality
    POLY_CHECK(key1 != key3);
    POLY_CHECK(!(key1 == key3));

    cout << "Operator tests passed!" << endl;
  }

  //----------------------------------------------------------------------------
  // Bit distribution test - ensure we use full 128-bit space
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Testing bit distribution ===" << endl;

    // Test that both lo and hi bits are used
    uint64_t lo_or = 0, hi_or = 0;
    for (unsigned i = 0; i < 1000; ++i) {
      uint64_t x = uint64_t(random01() * (1ULL << 41));
      uint64_t y = uint64_t(random01() * (1ULL << 41));
      uint64_t z = uint64_t(random01() * (1ULL << 41));

      auto key = HashKey3D::hash(Point3<uint64_t>(x, y, z));
      lo_or |= (uint64_t)key;
      hi_or |= (uint64_t)(key >> 64);
    }

    // Count bits set in lo and hi
    unsigned lo_bits = __builtin_popcountll(lo_or);
    unsigned hi_bits = __builtin_popcountll(hi_or);

    cout << "  Bits used: lo=" << lo_bits << "/64, hi=" << hi_bits << "/64" << endl;

    POLY_CHECK2(lo_bits > 50, "Too few bits used in lo word: " << lo_bits);
    POLY_CHECK2(hi_bits > 50, "Too few bits used in hi word: " << hi_bits);

    cout << "Bit distribution test passed!" << endl;
  }

  //----------------------------------------------------------------------------
  // Flag bit tests
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Testing flag bit functionality ===" << endl;

    // Test that freshly hashed points have flag disabled (inner box)
    auto key1 = HashKey3D::hash(Point3<uint64_t>(100, 200, 300));
    POLY_CHECK2(!HashKey3D::getOuterFlag(key1),
                "Freshly hashed point should have outer flag disabled");

    // Test enabling the flag
    HashKey3D::enableOuterFlag(key1);
    POLY_CHECK2(HashKey3D::getOuterFlag(key1),
                "Flag should be enabled after enableOuterFlag()");

    // Test disabling the flag
    HashKey3D::disableOuterFlag(key1);
    POLY_CHECK2(!HashKey3D::getOuterFlag(key1),
                "Flag should be disabled after disableOuterFlag()");

    // Test that flag doesn't affect unhashing
    Point3<uint64_t> p(12345, 67890, 11111);
    auto key2 = HashKey3D::hash(p);
    auto key3 = key2;
    HashKey3D::enableOuterFlag(key3);

    Point3<uint64_t> p2 = HashKey3D::unhash(key2);
    Point3<uint64_t> p3 = HashKey3D::unhash(key3);

    POLY_CHECK2(p2.x == p.x && p2.y == p.y && p2.z == p.z,
                "Unhashing with flag disabled should recover original point");
    POLY_CHECK2(p3.x == p.x && p3.y == p.y && p3.z == p.z,
                "Unhashing with flag enabled should recover original point");

    // Test that two hashes differ only in flag bit
    auto keyInner = HashKey3D::hash(Point3<uint64_t>(500, 750, 999));
    auto keyOuter = keyInner;
    HashKey3D::enableOuterFlag(keyOuter);

    POLY_CHECK2(keyInner != keyOuter,
                "Hashes with different flags should not be equal");
    POLY_CHECK2((keyInner ^ keyOuter) == HashKey3D::FlagMask(),
                "Hashes should differ only in the flag bit");

    // Test flag bit position
    POLY_CHECK2(HashKey3D::flagBit() == 127,
                "Flag bit should be at position 127 for 3D");

    // Test that flag mask has only one bit set
    unsigned __int128 mask = HashKey3D::FlagMask();
    unsigned bitCount = 0;
    unsigned __int128 temp = mask;
    while (temp) {
      bitCount += temp & 1;
      temp >>= 1;
    }
    POLY_CHECK2(bitCount == 1,
                "Flag mask should have exactly one bit set, found " << bitCount);

    // Test flag operations are idempotent
    auto key4 = HashKey3D::hash(Point3<uint64_t>(111, 222, 333));
    HashKey3D::enableOuterFlag(key4);
    HashKey3D::enableOuterFlag(key4);  // Enable twice
    POLY_CHECK2(HashKey3D::getOuterFlag(key4),
                "Double enable should still have flag set");

    HashKey3D::disableOuterFlag(key4);
    HashKey3D::disableOuterFlag(key4);  // Disable twice
    POLY_CHECK2(!HashKey3D::getOuterFlag(key4),
                "Double disable should still have flag clear");

    // Test that flag doesn't interfere with spatial locality
    // Hash the same point with and without flag - should sort adjacently
    Point3<uint64_t> testPt(1000, 2000, 3000);
    auto hashInner = HashKey3D::hash(testPt);
    auto hashOuter = hashInner;
    HashKey3D::enableOuterFlag(hashOuter);

    // XOR should only have the flag bit set
    unsigned __int128 diff = hashInner ^ hashOuter;
    POLY_CHECK2(diff == HashKey3D::FlagMask(),
                "Flag operations should only affect the flag bit");

    cout << "Flag bit tests passed!" << endl;
  }

  cout << "\n=== All HashKey3D tests passed! ===" << endl;

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
