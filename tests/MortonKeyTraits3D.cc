// Unit tests for HashKey3D Morton encoding/decoding (128-bit)
//
// NOTE: HashKey3D uses signed integers (int64_t) for coordinates,
// but only non-negative values (>= 0) are used in practice.

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <cassert>
#include <cstdint>
#include <algorithm>
#include <iomanip>

#include "polytope.hh"
#include "MortonKeyTraits.hh"
#include "Point.hh"

#include "Communicator.hh" 

using namespace std;
using namespace polytope;

namespace {
inline
double random01() {
  return double(rand())/RAND_MAX;
}

//------------------------------------------------------------------------------
// Type alias for HashKey3D
//------------------------------------------------------------------------------
using MortonKeyTraits3D = MortonKeyTraits<3>;
using Key = MortonKey<3>;
using IntType = QuantizedCoordinate<3>;
using IntPoint = QuantizedPoint<3>;

//------------------------------------------------------------------------------
// Test that encode/decode are inverses
//------------------------------------------------------------------------------
void testRoundTrip(const IntPoint& p, const string& label) {
  auto key = MortonKeyTraits3D::encode(p);
  auto p2 = MortonKeyTraits3D::decode(key);

  POLY_CHECK2(p.x == p2.x && p.y == p2.y && p.z == p2.z,
              label << ": Round-trip failed for (" << p.x << ", " << p.y << ", " << p.z
              << ") -> (" << p2.x << ", " << p2.y << ", " << p2.z << ")");
}

//------------------------------------------------------------------------------
// Test that identical points produce identical hashes
//------------------------------------------------------------------------------
void testIdenticalPoints(const IntPoint& p) {
  auto key1 = MortonKeyTraits3D::encode(p);
  auto key2 = MortonKeyTraits3D::encode(p);

  POLY_CHECK2(key1 == key2,
              "Identical points should hash identically: (" << p.x << ", " << p.y
              << ", " << p.z << ")");
}

//------------------------------------------------------------------------------
// Test that different points produce different hashes
//------------------------------------------------------------------------------
void testDistinctPoints(const IntPoint& p1, const IntPoint& p2) {
  if (p1.x == p2.x && p1.y == p2.y && p1.z == p2.z) return;

  auto key1 = MortonKeyTraits3D::encode(p1);
  auto key2 = MortonKeyTraits3D::encode(p2);

  POLY_CHECK2(key1 != key2,
              "Distinct points should hash differently: (" << p1.x << ", " << p1.y
              << ", " << p1.z << ") vs (" << p2.x << ", " << p2.y << ", " << p2.z << ")");
}

//------------------------------------------------------------------------------
// Test uniqueness - N distinct points produce N distinct hashes
//------------------------------------------------------------------------------
void testUniqueness(const vector<IntPoint>& points) {
  map<Key, vector<IntType>> collisions;

  for (unsigned i = 0; i < points.size(); ++i) {
    auto key = MortonKeyTraits3D::encode(points[i]);
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
  vector<IntPoint> points;
  points.reserve(gridSize * gridSize * gridSize);

  for (unsigned i = 0; i < gridSize; ++i) {
    for (unsigned j = 0; j < gridSize; ++j) {
      for (unsigned k = 0; k < gridSize; ++k) {
        points.push_back(IntPoint(i, j, k));
      }
    }
  }

  // Hash and sort by Morton code
  vector<pair<Key, unsigned>> hashed;
  for (unsigned i = 0; i < points.size(); ++i) {
    auto key = MortonKeyTraits3D::encode(points[i]);
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
  IntType zero = 0;
  testRoundTrip(IntPoint(zero, zero, zero), "Origin");

  // Max values for coordinate type (42 bits per dimension)
  const auto maxVal = MortonKeyTraits3D::maxCoordinate();
  testRoundTrip(IntPoint(maxVal, maxVal, maxVal), "Max coordinates");
  testRoundTrip(IntPoint(maxVal, zero, zero), "Max X");
  testRoundTrip(IntPoint(zero, maxVal, zero), "Max Y");
  testRoundTrip(IntPoint(zero, zero, maxVal), "Max Z");
  testRoundTrip(IntPoint(maxVal, maxVal, zero), "Max X,Y");
  testRoundTrip(IntPoint(maxVal, zero, maxVal), "Max X,Z");
  testRoundTrip(IntPoint(zero, maxVal, maxVal), "Max Y,Z");

  // Powers of 2
  for (auto bit = 0; bit < std::min(20, MortonKeyTraits3D::bitsPerCoordinate); ++bit) {
    auto val = IntType(1) << bit;
    testRoundTrip(IntPoint(val, zero, zero), "Power of 2 in X");
    testRoundTrip(IntPoint(zero, val, zero), "Power of 2 in Y");
    testRoundTrip(IntPoint(zero, zero, val), "Power of 2 in Z");
    testRoundTrip(IntPoint(val, val, val), "Power of 2 in all");
  }

  // Small integer grid
  for (auto x = 0; x < 8; ++x) {
    for (auto y = 0; y < 8; ++y) {
      for (auto z = 0; z < 8; ++z) {
        testRoundTrip(IntPoint(x, y, z), "Small integers");
        testIdenticalPoints(IntPoint(x, y, z));

        if (x + 1 < 8) {
          testDistinctPoints(IntPoint(x, y, z),
                             IntPoint(x+1, y, z));
        }
        if (y + 1 < 8) {
          testDistinctPoints(IntPoint(x, y, z),
                             IntPoint(x, y+1, z));
        }
        if (z + 1 < 8) {
          testDistinctPoints(IntPoint(x, y, z),
                             IntPoint(x, y, z+1));
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
  IntPoint p(0, 0, 0);

  // Set high bit in each dimension
  p.x = IntType(1) << 41;
  p.y = IntType(1) << 41;
  p.z = IntType(1) << 41;

  auto key = MortonKeyTraits3D::encode(p);
  auto p2 = MortonKeyTraits3D::decode(key);

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

  auto& comm = Communicator::instance();
  comm.init(argc, argv);

  //----------------------------------------------------------------------------
  // Test 128-bit structure
  //----------------------------------------------------------------------------
  test128BitStructure();

  //----------------------------------------------------------------------------
  // Test with 42-bit range (non-negative values)
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Testing HashKey3D with 42-bit range ===" << endl;

    // Edge cases
    testEdgeCases();

    // Random points in non-negative 42-bit range
    cout << "Testing random points (42-bit non-negative range)..." << endl;
    const unsigned n = 800;
    vector<IntPoint> points;
    const auto range = MortonKeyTraits3D::maxCoordinate();
    for (unsigned i = 0; i < n; ++i) {
      auto x = IntType(random01() * range);
      auto y = IntType(random01() * range);
      auto z = IntType(random01() * range);
      points.push_back(IntPoint(x, y, z));
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
    vector<IntPoint> points;
    for (unsigned i = 0; i < n; ++i) {
      auto x = IntType(random01() * (IntType(1) << 20));
      auto y = IntType(random01() * (IntType(1) << 20));
      auto z = IntType(random01() * (IntType(1) << 20));
      points.push_back(IntPoint(x, y, z));
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

    const auto range = MortonKeyTraits3D::maxCoordinate();
    for (unsigned i = 0; i < nOps; ++i) {
      auto x = IntType(random01() * range);
      auto y = IntType(random01() * range);
      auto z = IntType(random01() * range);

      auto key = MortonKeyTraits3D::encode(IntPoint(x, y, z));
      auto recovered = MortonKeyTraits3D::decode(key);

      POLY_CHECK(recovered.x == x && recovered.y == y && recovered.z == z);
    }
    cout << "Stress test passed!" << endl;
  }

  //----------------------------------------------------------------------------
  // Operator tests
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Testing operators ===" << endl;

    auto key1 = MortonKeyTraits3D::encode(IntPoint(42, 73, 101));
    auto key2 = MortonKeyTraits3D::encode(IntPoint(42, 73, 101));
    auto key3 = MortonKeyTraits3D::encode(IntPoint(101, 73, 42));

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
    const auto range = MortonKeyTraits3D::maxCoordinate();
    for (unsigned i = 0; i < 1000; ++i) {
      auto x = IntType(random01() * range);
      auto y = IntType(random01() * range);
      auto z = IntType(random01() * range);

      auto key = MortonKeyTraits3D::encode(IntPoint(x, y, z));
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

  cout << "\n=== All HashKey3D tests passed! ===" << endl;

  comm.finalize();
  return 0;
}
