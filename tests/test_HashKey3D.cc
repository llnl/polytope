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
// Helper to compare 128-bit hash keys
//------------------------------------------------------------------------------
bool hashEqual(const HashKey3D& a, const HashKey3D& b) {
  return a.lo() == b.lo() && a.hi() == b.hi();
}

//------------------------------------------------------------------------------
// Test that interleave/deinterleave are inverses
//------------------------------------------------------------------------------
template<typename CoordType>
void testRoundTrip(const Point3<CoordType>& p, const string& label) {
  HashKey3D key;
  key.interleave(p);
  Point3<CoordType> p2 = key.deinterleave<CoordType>();

  POLY_CHECK2(p.x == p2.x && p.y == p2.y && p.z == p2.z,
              label << ": Round-trip failed for (" << p.x << ", " << p.y << ", " << p.z
              << ") -> (" << p2.x << ", " << p2.y << ", " << p2.z << ")");
}

//------------------------------------------------------------------------------
// Test that identical points produce identical hashes
//------------------------------------------------------------------------------
template<typename CoordType>
void testIdenticalPoints(const Point3<CoordType>& p) {
  HashKey3D key1, key2;
  key1.interleave(p);
  key2.interleave(p);

  POLY_CHECK2(key1 == key2,
              "Identical points should hash identically: (" << p.x << ", " << p.y
              << ", " << p.z << ")");
  POLY_CHECK2(hashEqual(key1, key2),
              "Hash values should match for identical points");
}

//------------------------------------------------------------------------------
// Test that different points produce different hashes
//------------------------------------------------------------------------------
template<typename CoordType>
void testDistinctPoints(const Point3<CoordType>& p1, const Point3<CoordType>& p2) {
  if (p1.x == p2.x && p1.y == p2.y && p1.z == p2.z) return;

  HashKey3D key1, key2;
  key1.interleave(p1);
  key2.interleave(p2);

  POLY_CHECK2(key1 != key2,
              "Distinct points should hash differently: (" << p1.x << ", " << p1.y
              << ", " << p1.z << ") vs (" << p2.x << ", " << p2.y << ", " << p2.z << ")");
}

//------------------------------------------------------------------------------
// Test uniqueness - N distinct points produce N distinct hashes
//------------------------------------------------------------------------------
template<typename CoordType>
void testUniqueness(const vector<Point3<CoordType>>& points) {
  map<pair<uint64_t, uint64_t>, vector<unsigned>> collisions;

  for (unsigned i = 0; i < points.size(); ++i) {
    HashKey3D key;
    key.interleave(points[i]);
    pair<uint64_t, uint64_t> h(key.lo(), key.hi());
    collisions[h].push_back(i);
  }

  // Report any collisions
  unsigned collisionCount = 0;
  for (const auto& kv : collisions) {
    if (kv.second.size() > 1) {
      collisionCount++;
      cout << "  Collision at hash (" << hex << kv.first.first << ", "
           << kv.first.second << dec << "): ";
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
template<typename CoordType>
void testSpatialLocality() {
  cout << "Testing spatial locality with 3D grid..." << endl;

  // Create a regular 3D grid of points (16x16x16 = 4096 points)
  const unsigned gridSize = 16;
  vector<Point3<CoordType>> points;
  points.reserve(gridSize * gridSize * gridSize);

  for (unsigned i = 0; i < gridSize; ++i) {
    for (unsigned j = 0; j < gridSize; ++j) {
      for (unsigned k = 0; k < gridSize; ++k) {
        points.push_back(Point3<CoordType>(CoordType(i), CoordType(j), CoordType(k)));
      }
    }
  }

  // Hash and sort by Morton code
  vector<pair<pair<uint64_t, uint64_t>, unsigned>> hashed;
  for (unsigned i = 0; i < points.size(); ++i) {
    HashKey3D key;
    key.interleave(points[i]);
    hashed.push_back(make_pair(make_pair(key.lo(), key.hi()), i));
  }

  // Sort by hash value (hi, then lo)
  sort(hashed.begin(), hashed.end(),
       [](const pair<pair<uint64_t, uint64_t>, unsigned>& a,
          const pair<pair<uint64_t, uint64_t>, unsigned>& b) {
         if (a.first.second != b.first.second)
           return a.first.second < b.first.second;
         return a.first.first < b.first.first;
       });

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
template<typename CoordType>
void testEdgeCases() {
  cout << "Testing edge cases..." << endl;

  // Origin
  testRoundTrip(Point3<CoordType>(0, 0, 0), "Origin");

  // Max values for coordinate type (42 bits per dimension)
  const CoordType maxVal = (CoordType(1) << (HashKey3D::bitsPerDim() - 1)) - 1;
  testRoundTrip(Point3<CoordType>(maxVal, maxVal, maxVal), "Max coordinates");
  testRoundTrip(Point3<CoordType>(maxVal, 0, 0), "Max X");
  testRoundTrip(Point3<CoordType>(0, maxVal, 0), "Max Y");
  testRoundTrip(Point3<CoordType>(0, 0, maxVal), "Max Z");
  testRoundTrip(Point3<CoordType>(maxVal, maxVal, 0), "Max X,Y");
  testRoundTrip(Point3<CoordType>(maxVal, 0, maxVal), "Max X,Z");
  testRoundTrip(Point3<CoordType>(0, maxVal, maxVal), "Max Y,Z");

  // Powers of 2
  for (unsigned bit = 0; bit < min(20u, HashKey3D::bitsPerDim()); ++bit) {
    CoordType val = CoordType(1) << bit;
    testRoundTrip(Point3<CoordType>(val, 0, 0), "Power of 2 in X");
    testRoundTrip(Point3<CoordType>(0, val, 0), "Power of 2 in Y");
    testRoundTrip(Point3<CoordType>(0, 0, val), "Power of 2 in Z");
    testRoundTrip(Point3<CoordType>(val, val, val), "Power of 2 in all");
  }

  // Small integer grid
  for (CoordType x = 0; x < 8; ++x) {
    for (CoordType y = 0; y < 8; ++y) {
      for (CoordType z = 0; z < 8; ++z) {
        testRoundTrip(Point3<CoordType>(x, y, z), "Small integers");
        testIdenticalPoints(Point3<CoordType>(x, y, z));

        if (x + 1 < 8) {
          testDistinctPoints(Point3<CoordType>(x, y, z),
                           Point3<CoordType>(x+1, y, z));
        }
        if (y + 1 < 8) {
          testDistinctPoints(Point3<CoordType>(x, y, z),
                           Point3<CoordType>(x, y+1, z));
        }
        if (z + 1 < 8) {
          testDistinctPoints(Point3<CoordType>(x, y, z),
                           Point3<CoordType>(x, y, z+1));
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
  HashKey3D key;
  Point3<uint64_t> p(0, 0, 0);

  // Set high bit in each dimension
  p.x = uint64_t(1) << 41;
  p.y = uint64_t(1) << 41;
  p.z = uint64_t(1) << 41;

  key.interleave(p);
  Point3<uint64_t> p2 = key.deinterleave<uint64_t>();

  POLY_CHECK2(p.x == p2.x && p.y == p2.y && p.z == p2.z,
              "Failed to represent 42-bit values: (" << p.x << ", " << p.y
              << ", " << p.z << ") -> (" << p2.x << ", " << p2.y << ", " << p2.z << ")");

  // Test that lo and hi are used
  cout << "  Hash of (2^41, 2^41, 2^41): lo=" << hex << key.lo()
       << " hi=" << key.hi() << dec << endl;
  POLY_CHECK2(key.hi() != 0, "High 64 bits should be used for 42-bit coordinates");

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
    testEdgeCases<uint64_t>();

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
    testSpatialLocality<uint64_t>();
  }

  //----------------------------------------------------------------------------
  // Test 32-bit unsigned integers
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Testing HashKey3D with uint32_t ===" << endl;

    const unsigned n = 1000;
    vector<Point3<uint32_t>> points;
    for (unsigned i = 0; i < n; ++i) {
      uint32_t x = uint32_t(random01() * (1ULL << 30));
      uint32_t y = uint32_t(random01() * (1ULL << 30));
      uint32_t z = uint32_t(random01() * (1ULL << 30));
      points.push_back(Point3<uint32_t>(x, y, z));
      testRoundTrip(points.back(), "Random 32-bit point");
    }

    testUniqueness(points);
    testSpatialLocality<uint32_t>();
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
    testSpatialLocality<uint64_t>();
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

      HashKey3D key;
      key.interleave(Point3<uint64_t>(x, y, z));
      Point3<uint64_t> recovered = key.deinterleave<uint64_t>();

      POLY_CHECK(recovered.x == x && recovered.y == y && recovered.z == z);
    }
    cout << "Stress test passed!" << endl;
  }

  //----------------------------------------------------------------------------
  // Operator tests
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Testing operators ===" << endl;

    HashKey3D key1, key2, key3;
    key1.interleave(Point3<uint64_t>(42, 73, 101));
    key2.interleave(Point3<uint64_t>(42, 73, 101));
    key3.interleave(Point3<uint64_t>(101, 73, 42));

    // Equality
    POLY_CHECK(key1 == key2);
    POLY_CHECK(!(key1 != key2));

    // Inequality
    POLY_CHECK(key1 != key3);
    POLY_CHECK(!(key1 == key3));

    // lo/hi accessors
    POLY_CHECK(key1.lo() == key2.lo());
    POLY_CHECK(key1.hi() == key2.hi());
    POLY_CHECK(key1.lo() != key3.lo() || key1.hi() != key3.hi());

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

      HashKey3D key;
      key.interleave(Point3<uint64_t>(x, y, z));
      lo_or |= key.lo();
      hi_or |= key.hi();
    }

    // Count bits set in lo and hi
    unsigned lo_bits = __builtin_popcountll(lo_or);
    unsigned hi_bits = __builtin_popcountll(hi_or);

    cout << "  Bits used: lo=" << lo_bits << "/64, hi=" << hi_bits << "/64" << endl;

    POLY_CHECK2(lo_bits > 50, "Too few bits used in lo word: " << lo_bits);
    POLY_CHECK2(hi_bits > 50, "Too few bits used in hi word: " << hi_bits);

    cout << "Bit distribution test passed!" << endl;
  }

  cout << "\n=== All HashKey3D tests passed! ===" << endl;

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
