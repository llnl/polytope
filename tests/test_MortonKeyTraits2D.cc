// Unit tests for HashKey2D Morton encoding/decoding
//
// NOTE: HashKey2D uses signed integers (int) for coordinates,
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
// Type aliases for 2D Morton encoding
//------------------------------------------------------------------------------
using MortonKeyTraits2D = MortonKeyTraits<2>;
using Key = MortonKey<2>;
using IntType = QuantizedCoordinate<2>;
using IntPoint = QuantizedPoint<2>;

//------------------------------------------------------------------------------
// Test that encode/decode are inverses
//------------------------------------------------------------------------------
void testRoundTrip(const IntPoint& p, const string& label) {
  auto key = MortonKeyTraits2D::encode(p);
  auto p2 = MortonKeyTraits2D::decode(key);

  POLY_CHECK2(p.x == p2.x && p.y == p2.y,
              label << ": Round-trip failed for (" << p.x << ", " << p.y << ") -> ("
              << p2.x << ", " << p2.y << ")");
}

//------------------------------------------------------------------------------
// Test that identical points produce identical hashes
//------------------------------------------------------------------------------
void testIdenticalPoints(const IntPoint& p) {
  auto key1 = MortonKeyTraits2D::encode(p);
  auto key2 = MortonKeyTraits2D::encode(p);

  POLY_CHECK2(key1 == key2,
              "Identical points should hash identically: (" << p.x << ", " << p.y << ")");
}

//------------------------------------------------------------------------------
// Test that different points produce different hashes
//------------------------------------------------------------------------------
void testDistinctPoints(const IntPoint& p1, const IntPoint& p2) {
  if (p1.x == p2.x && p1.y == p2.y) return;

  auto key1 = MortonKeyTraits2D::encode(p1);
  auto key2 = MortonKeyTraits2D::encode(p2);

  POLY_CHECK2(key1 != key2,
              "Distinct points should hash differently: (" << p1.x << ", " << p1.y
              << ") vs (" << p2.x << ", " << p2.y << ")");
}

//------------------------------------------------------------------------------
// Test spatial locality - Morton curves place nearby points near each other
// Uses a grid of points to ensure meaningful locality testing
//------------------------------------------------------------------------------
void testSpatialLocality() {
  cout << "Testing spatial locality with grid..." << endl;

  // Create a regular grid of points (32x32 = 1024 points)
  const unsigned gridSize = 32;
  vector<IntPoint> points;
  points.reserve(gridSize * gridSize);

  for (unsigned i = 0; i < gridSize; ++i) {
    for (unsigned j = 0; j < gridSize; ++j) {
      points.push_back(IntPoint(i, j));
    }
  }

  // Hash and sort by Morton code
  vector<pair<Key, unsigned>> hashed;
  for (unsigned i = 0; i < points.size(); ++i) {
    auto key = MortonKeyTraits2D::encode(points[i]);
    hashed.push_back(make_pair(key, i));
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
void testUniqueness(const vector<IntPoint>& points) {
  set<Key> hashes;
  map<Key, vector<IntType>> collisions;

  for (unsigned i = 0; i < points.size(); ++i) {
    auto key = MortonKeyTraits2D::encode(points[i]);

    if (hashes.count(key)) {
      collisions[key].push_back(i);
    } else {
      hashes.insert(key);
      collisions[key] = {IntType(i)};
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
void testEdgeCases() {
  cout << "Testing edge cases..." << endl;

  // Origin
  IntType zero = 0;
  testRoundTrip(IntPoint(zero, zero), "Origin");

  // Max values for the coordinate type
  const auto maxVal = MortonKeyTraits2D::maxCoordinate();
  testRoundTrip(IntPoint(maxVal, maxVal), "Max coordinates");
  testRoundTrip(IntPoint(maxVal, zero), "Max X, zero Y");
  testRoundTrip(IntPoint(zero, maxVal), "Zero X, max Y");

  // Powers of 2
  for (auto bit = 0; bit < std::min(20, MortonKeyTraits2D::bitsPerCoordinate); ++bit) {
    auto val = IntType(1) << bit;
    testRoundTrip(IntPoint(val, zero), "Power of 2 in X");
    testRoundTrip(IntPoint(zero, val), "Power of 2 in Y");
    testRoundTrip(IntPoint(val, val), "Power of 2 in both");
  }

  // Small coordinates
  for (auto x = 0; x < 10; ++x) {
    for (auto y = 0; y < 10; ++y) {
      testRoundTrip(IntPoint(x, y), "Small integers");
      testIdenticalPoints(IntPoint(x, y));

      if (x + 1 < 10) {
        testDistinctPoints(IntPoint(x, y), IntPoint(x+1, y));
      }
      if (y + 1 < 10) {
        testDistinctPoints(IntPoint(x, y), IntPoint(x, y+1));
      }
    }
  }
}

} // anonymous namespace

//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

  auto& comm = Communicator::instance();
  comm.init(argc, argv);

  //----------------------------------------------------------------------------
  // Test with 30-bit range (non-negative values)
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Testing HashKey2D with 30-bit range ===" << endl;

    // Edge cases
    testEdgeCases();

    // Random points in non-negative range
    cout << "Testing random points (30-bit non-negative range)..." << endl;
    const unsigned n = 1000;
    vector<IntPoint> points;
    const auto range = MortonKeyTraits2D::maxCoordinate();
    for (unsigned i = 0; i < n; ++i) {
      auto x = IntType(random01() * range);
      auto y = IntType(random01() * range);
      points.push_back(IntPoint(x, y));
      testRoundTrip(points.back(), "Random point");
    }

    // Uniqueness test
    testUniqueness(points);

    // Spatial locality test (uses its own grid)
    testSpatialLocality();
  }

  //----------------------------------------------------------------------------
  // Test with smaller range (16-bit values)
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Testing HashKey2D with 16-bit range ===" << endl;

    const unsigned n = 1000;
    vector<IntPoint> points;
    for (unsigned i = 0; i < n; ++i) {
      auto x = IntType(random01() * 65536);
      auto y = IntType(random01() * 65536);
      points.push_back(IntPoint(x, y));
      testRoundTrip(points.back(), "Random 16-bit point");
    }

    testUniqueness(points);
    testSpatialLocality();
  }

  //----------------------------------------------------------------------------
  // Stress test - Ensure no crashes with many operations
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Stress test ===" << endl;
    const unsigned nOps = 10000;
    cout << "Performing " << nOps << " hash/unhash operations..." << endl;

    const auto range = MortonKeyTraits2D::maxCoordinate();
    for (unsigned i = 0; i < nOps; ++i) {
      auto x = IntType(random01() * range);
      auto y = IntType(random01() * range);

      auto key = MortonKeyTraits2D::encode(IntPoint(x, y));
      auto recovered = MortonKeyTraits2D::decode(key);

      POLY_CHECK(recovered.x == x && recovered.y == y);
    }
    cout << "Stress test passed!" << endl;
  }

  //----------------------------------------------------------------------------
  // Operator tests
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Testing operators ===" << endl;

    auto key1 = MortonKeyTraits2D::encode(IntPoint(42, 73));
    auto key2 = MortonKeyTraits2D::encode(IntPoint(42, 73));
    auto key3 = MortonKeyTraits2D::encode(IntPoint(73, 42));

    // Equality
    POLY_CHECK(key1 == key2);
    POLY_CHECK(!(key1 != key2));

    // Inequality
    POLY_CHECK(key1 != key3);
    POLY_CHECK(!(key1 == key3));

    cout << "Operator tests passed!" << endl;
  }

  cout << "\n=== All HashKey2D tests passed! ===" << endl;

  comm.finalize();
  return 0;
}
