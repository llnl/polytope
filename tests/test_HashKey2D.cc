// Unit tests for HashKey2D Morton encoding/decoding
//
// NOTE: HashKey2D only supports non-negative coordinates.
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
// Test that interleave/deinterleave are inverses
//------------------------------------------------------------------------------
template<typename CoordType>
void testRoundTrip(const Point2<CoordType>& p, const string& label) {
  HashKey2D key;
  key.interleave(p);
  Point2<CoordType> p2 = key.deinterleave<CoordType>();

  POLY_CHECK2(p.x == p2.x && p.y == p2.y,
              label << ": Round-trip failed for (" << p.x << ", " << p.y << ") -> ("
              << p2.x << ", " << p2.y << ")");
}

//------------------------------------------------------------------------------
// Test that identical points produce identical hashes
//------------------------------------------------------------------------------
template<typename CoordType>
void testIdenticalPoints(const Point2<CoordType>& p) {
  HashKey2D key1, key2;
  key1.interleave(p);
  key2.interleave(p);

  POLY_CHECK2(key1 == key2,
              "Identical points should hash identically: (" << p.x << ", " << p.y << ")");
  POLY_CHECK2(key1.value() == key2.value(),
              "Hash values should match for identical points");
}

//------------------------------------------------------------------------------
// Test that different points produce different hashes
//------------------------------------------------------------------------------
template<typename CoordType>
void testDistinctPoints(const Point2<CoordType>& p1, const Point2<CoordType>& p2) {
  if (p1.x == p2.x && p1.y == p2.y) return;

  HashKey2D key1, key2;
  key1.interleave(p1);
  key2.interleave(p2);

  POLY_CHECK2(key1 != key2,
              "Distinct points should hash differently: (" << p1.x << ", " << p1.y
              << ") vs (" << p2.x << ", " << p2.y << ")");
}

//------------------------------------------------------------------------------
// Test spatial locality - Morton curves place nearby points near each other
// Uses a grid of points to ensure meaningful locality testing
//------------------------------------------------------------------------------
template<typename CoordType>
void testSpatialLocality() {
  cout << "Testing spatial locality with grid..." << endl;

  // Create a regular grid of points (32x32 = 1024 points)
  const unsigned gridSize = 32;
  vector<Point2<CoordType>> points;
  points.reserve(gridSize * gridSize);

  for (unsigned i = 0; i < gridSize; ++i) {
    for (unsigned j = 0; j < gridSize; ++j) {
      points.push_back(Point2<CoordType>(CoordType(i), CoordType(j)));
    }
  }

  // Hash and sort by Morton code
  vector<pair<uint64_t, unsigned>> hashed;
  for (unsigned i = 0; i < points.size(); ++i) {
    HashKey2D key;
    key.interleave(points[i]);
    hashed.push_back(make_pair(key.value(), i));
  }
  sort(hashed.begin(), hashed.end());

  // Compute average distance between consecutive points in Morton order
  double totalDist = 0.0;
  unsigned nearbyCount = 0;
  for (unsigned i = 0; i + 1 < hashed.size(); ++i) {
    unsigned idx1 = hashed[i].second;
    unsigned idx2 = hashed[i+1].second;

    double dx = double(points[idx1].x > points[idx2].x ?
                       points[idx1].x - points[idx2].x : points[idx2].x - points[idx1].x);
    double dy = double(points[idx1].y > points[idx2].y ?
                       points[idx1].y - points[idx2].y : points[idx2].y - points[idx1].y);
    double dist = dx + dy; // Manhattan distance
    totalDist += dist;

    // Points within small distance are "nearby"
    if (dist <= 4.0) nearbyCount++;
  }

  double avgDist = totalDist / (hashed.size() - 1);
  double localityRatio = double(nearbyCount) / (hashed.size() - 1);

  cout << "  Average Manhattan distance between consecutive hashes: " << avgDist << endl;
  cout << "  Spatial locality ratio (dist <= 4): " << localityRatio << " ("
       << nearbyCount << "/" << (hashed.size() - 1) << ")" << endl;

  // Morton curves should have good locality - average distance should be small
  // For a 32x32 grid, random pairing would give avg distance ~21
  POLY_CHECK2(avgDist < 10.0,
              "Average distance too high: " << avgDist << " (Morton locality poor)");
  POLY_CHECK2(localityRatio > 0.4,
              "Spatial locality ratio too low: " << localityRatio);
}

//------------------------------------------------------------------------------
// Test uniqueness - N distinct points produce N distinct hashes
//------------------------------------------------------------------------------
template<typename CoordType>
void testUniqueness(const vector<Point2<CoordType>>& points) {
  set<uint64_t> hashes;
  map<uint64_t, vector<unsigned>> collisions;

  for (unsigned i = 0; i < points.size(); ++i) {
    HashKey2D key;
    key.interleave(points[i]);
    uint64_t h = key.value();

    if (hashes.count(h)) {
      collisions[h].push_back(i);
    } else {
      hashes.insert(h);
      collisions[h] = {i};
    }
  }

  // Report any collisions
  unsigned collisionCount = 0;
  for (const auto& kv : collisions) {
    if (kv.second.size() > 1) {
      collisionCount++;
      cout << "  Collision at hash " << hex << kv.first << dec << ": ";
      for (unsigned idx : kv.second) {
        cout << "(" << points[idx].x << "," << points[idx].y << ") ";
      }
      cout << endl;
    }
  }

  POLY_CHECK2(collisionCount == 0,
              "Found " << collisionCount << " hash collisions among "
              << points.size() << " distinct points");
}

//------------------------------------------------------------------------------
// Test edge cases
//------------------------------------------------------------------------------
template<typename CoordType>
void testEdgeCases() {
  cout << "Testing edge cases..." << endl;

  // Origin
  testRoundTrip(Point2<CoordType>(0, 0), "Origin");

  // Max values for the coordinate type
  const CoordType maxVal = (CoordType(1) << (HashKey2D::bitsPerDim() - 1)) - 1;
  testRoundTrip(Point2<CoordType>(maxVal, maxVal), "Max coordinates");
  testRoundTrip(Point2<CoordType>(maxVal, 0), "Max X, zero Y");
  testRoundTrip(Point2<CoordType>(0, maxVal), "Zero X, max Y");

  // Powers of 2
  for (unsigned bit = 0; bit < min(20u, HashKey2D::bitsPerDim()); ++bit) {
    CoordType val = CoordType(1) << bit;
    testRoundTrip(Point2<CoordType>(val, 0), "Power of 2 in X");
    testRoundTrip(Point2<CoordType>(0, val), "Power of 2 in Y");
    testRoundTrip(Point2<CoordType>(val, val), "Power of 2 in both");
  }

  // Adjacent coordinates
  for (CoordType x = 0; x < 10; ++x) {
    for (CoordType y = 0; y < 10; ++y) {
      testRoundTrip(Point2<CoordType>(x, y), "Small integers");
      testIdenticalPoints(Point2<CoordType>(x, y));

      if (x + 1 < 10) {
        testDistinctPoints(Point2<CoordType>(x, y), Point2<CoordType>(x+1, y));
      }
      if (y + 1 < 10) {
        testDistinctPoints(Point2<CoordType>(x, y), Point2<CoordType>(x, y+1));
      }
    }
  }
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
  // Test 32-bit unsigned integers
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Testing HashKey2D with uint32_t ===" << endl;

    // Edge cases
    testEdgeCases<uint32_t>();

    // Random points in full range
    cout << "Testing random points (full range)..." << endl;
    const unsigned n = 1000;
    vector<Point2<uint32_t>> points;
    for (unsigned i = 0; i < n; ++i) {
      uint32_t x = uint32_t(random01() * (uint64_t(1) << 31));
      uint32_t y = uint32_t(random01() * (uint64_t(1) << 31));
      points.push_back(Point2<uint32_t>(x, y));
      testRoundTrip(points.back(), "Random point");
    }

    // Uniqueness test
    testUniqueness(points);

    // Spatial locality test (uses its own grid)
    testSpatialLocality<uint32_t>();
  }

  //----------------------------------------------------------------------------
  // Test with smaller range (16-bit values in uint32_t)
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Testing HashKey2D with 16-bit range ===" << endl;

    const unsigned n = 1000;
    vector<Point2<uint32_t>> points;
    for (unsigned i = 0; i < n; ++i) {
      uint32_t x = uint32_t(random01() * 65536);
      uint32_t y = uint32_t(random01() * 65536);
      points.push_back(Point2<uint32_t>(x, y));
      testRoundTrip(points.back(), "Random 16-bit point");
    }

    testUniqueness(points);
    testSpatialLocality<uint32_t>();
  }

  //----------------------------------------------------------------------------
  // Test 64-bit unsigned integers (will use lower 32 bits)
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Testing HashKey2D with uint64_t (32-bit values) ===" << endl;

    const unsigned n = 500;
    vector<Point2<uint64_t>> points;
    for (unsigned i = 0; i < n; ++i) {
      // Keep in 32-bit range
      uint64_t x = uint64_t(random01() * (uint64_t(1) << 30));
      uint64_t y = uint64_t(random01() * (uint64_t(1) << 30));
      points.push_back(Point2<uint64_t>(x, y));
      testRoundTrip(points.back(), "Random 64-bit point");
    }

    testUniqueness(points);
    testSpatialLocality<uint64_t>();
  }

  //----------------------------------------------------------------------------
  // Stress test - Ensure no crashes with many operations
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Stress test ===" << endl;
    const unsigned nOps = 10000;
    cout << "Performing " << nOps << " hash/unhash operations..." << endl;

    for (unsigned i = 0; i < nOps; ++i) {
      uint32_t x = uint32_t(random01() * UINT32_MAX);
      uint32_t y = uint32_t(random01() * UINT32_MAX);

      HashKey2D key;
      key.interleave(Point2<uint32_t>(x, y));
      Point2<uint32_t> recovered = key.deinterleave<uint32_t>();

      POLY_CHECK(recovered.x == x && recovered.y == y);
    }
    cout << "Stress test passed!" << endl;
  }

  //----------------------------------------------------------------------------
  // Operator tests
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Testing operators ===" << endl;

    HashKey2D key1, key2, key3;
    key1.interleave(Point2<uint32_t>(42, 73));
    key2.interleave(Point2<uint32_t>(42, 73));
    key3.interleave(Point2<uint32_t>(73, 42));

    // Equality
    POLY_CHECK(key1 == key2);
    POLY_CHECK(!(key1 != key2));

    // Inequality
    POLY_CHECK(key1 != key3);
    POLY_CHECK(!(key1 == key3));

    // Value accessor
    POLY_CHECK(key1.value() == key2.value());
    POLY_CHECK(key1.value() != key3.value());

    cout << "Operator tests passed!" << endl;
  }

  cout << "\n=== All HashKey2D tests passed! ===" << endl;

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
