// Unit tests for 2D Morton encoding/decoding
//
// Only non-negative coordinates are encoded in practice.

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <cassert>
#include <cstdint>
#include <algorithm>
#include <iomanip>
#include <sstream>

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
using Traits = MortonKeyTraits<2>;
using Key = MortonKey<2>;
using Coordinate = QuantizedCoordinate<2>;
using PointType = QuantizedPoint<2>;

std::string keyString(const Key key) {
#ifdef POLYTOPE_ENABLE_HIBIT2D
  const auto lower = static_cast<std::uint64_t>(key);
  const auto upper = static_cast<std::uint64_t>(key >> 64);
  std::ostringstream result;
  result << std::hex << upper << ":" << lower;
  return result.str();
#else
  std::ostringstream result;
  result << std::hex << key;
  return result.str();
#endif
}

//------------------------------------------------------------------------------
// Test that encode/decode are inverses
//------------------------------------------------------------------------------
void testRoundTrip(const PointType& p, const string& label) {
  auto key = Traits::encode(p);
  auto p2 = Traits::decode(key);

  POLY_CHECK2(p.x == p2.x && p.y == p2.y,
              label << ": Round-trip failed for (" << p.x << ", " << p.y << ") -> ("
              << p2.x << ", " << p2.y << ")");
}

//------------------------------------------------------------------------------
// Test that identical points produce identical hashes
//------------------------------------------------------------------------------
void testIdenticalPoints(const PointType& p) {
  auto key1 = Traits::encode(p);
  auto key2 = Traits::encode(p);

  POLY_CHECK2(key1 == key2,
              "Identical points should hash identically: (" << p.x << ", " << p.y << ")");
}

//------------------------------------------------------------------------------
// Test that different points produce different hashes
//------------------------------------------------------------------------------
void testDistinctPoints(const PointType& p1, const PointType& p2) {
  if (p1.x == p2.x && p1.y == p2.y) return;

  auto key1 = Traits::encode(p1);
  auto key2 = Traits::encode(p2);

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
  vector<PointType> points;
  points.reserve(gridSize * gridSize);

  for (unsigned i = 0; i < gridSize; ++i) {
    for (unsigned j = 0; j < gridSize; ++j) {
      points.push_back(PointType(i, j));
    }
  }

  // Hash and sort by Morton code
  vector<pair<Key, unsigned>> hashed;
  for (unsigned i = 0; i < points.size(); ++i) {
    auto key = Traits::encode(points[i]);
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
void testUniqueness(const vector<PointType>& points) {
  set<Key> hashes;
  map<Key, vector<Coordinate>> collisions;

  for (unsigned i = 0; i < points.size(); ++i) {
    auto key = Traits::encode(points[i]);

    if (hashes.count(key)) {
      collisions[key].push_back(i);
    } else {
      hashes.insert(key);
      collisions[key] = {Coordinate(i)};
    }
  }

  // Report any collisions
  unsigned collisionCount = 0;
  for (const auto& kv : collisions) {
    if (kv.second.size() > 1) {
      collisionCount++;
      cout << "  Collision at key " << keyString(kv.first) << ": ";
      for (unsigned idx : kv.second) {
        cout << "(" << points[idx].x << "," << points[idx].y << ") ";
      }
      cout << endl;
    }
  }

  POLY_CHECK2(collisionCount == 0,
              "Found " << collisionCount << " key collisions among "
              << points.size() << " distinct points");
}

//------------------------------------------------------------------------------
// Test edge cases
//------------------------------------------------------------------------------
void testEdgeCases() {
  cout << "Testing edge cases..." << endl;

  // Origin
  Coordinate zero = 0;
  testRoundTrip(PointType(zero, zero), "Origin");

  // Max values for the coordinate type
  const auto maxVal = Traits::maxCoordinate();
  testRoundTrip(PointType(maxVal, maxVal), "Max coordinates");
  testRoundTrip(PointType(maxVal, zero), "Max X, zero Y");
  testRoundTrip(PointType(zero, maxVal), "Zero X, max Y");

  // Powers of 2
  for (auto bit = 0; bit < Traits::bitsPerCoordinate; ++bit) {
    const auto value = Coordinate{1} << bit;
    testRoundTrip(PointType(value, zero), "Power of 2 in X");
    testRoundTrip(PointType(zero, value), "Power of 2 in Y");
    testRoundTrip(PointType(value, value), "Power of 2 in both");
  }

  // Small coordinates
  for (auto x = 0; x < 10; ++x) {
    for (auto y = 0; y < 10; ++y) {
      testRoundTrip(PointType(x, y), "Small integers");
      testIdenticalPoints(PointType(x, y));

      if (x + 1 < 10) {
        testDistinctPoints(PointType(x, y), PointType(x+1, y));
      }
      if (y + 1 < 10) {
        testDistinctPoints(PointType(x, y), PointType(x, y+1));
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
  // Test the full configured non-negative coordinate range.
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Testing 2D Morton keys with "
         << Traits::bitsPerCoordinate << " coordinate bits ===" << endl;

    // Edge cases
    testEdgeCases();

    // Random points in non-negative range
    cout << "Testing random points (30-bit non-negative range)..." << endl;
    const unsigned n = 1000;
    vector<PointType> points;
    const auto range = Traits::maxCoordinate();
    for (unsigned i = 0; i < n; ++i) {
      auto x = Coordinate(random01() * range);
      auto y = Coordinate(random01() * range);
      points.push_back(PointType(x, y));
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
    cout << "\n=== Testing 2D Morton keys with 16-bit values ===" << endl;

    const unsigned n = 1000;
    vector<PointType> points;
    for (unsigned i = 0; i < n; ++i) {
      auto x = Coordinate(random01() * 65536);
      auto y = Coordinate(random01() * 65536);
      points.push_back(PointType(x, y));
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
    cout << "Performing " << nOps << " encode/decode operations..." << endl;

    const auto range = Traits::maxCoordinate();
    for (unsigned i = 0; i < nOps; ++i) {
      auto x = Coordinate(random01() * range);
      auto y = Coordinate(random01() * range);

      auto key = Traits::encode(PointType(x, y));
      auto recovered = Traits::decode(key);

      POLY_CHECK(recovered.x == x && recovered.y == y);
    }
    cout << "Stress test passed!" << endl;
  }

  //----------------------------------------------------------------------------
  // Operator tests
  //----------------------------------------------------------------------------
  {
    cout << "\n=== Testing operators ===" << endl;

    auto key1 = Traits::encode(PointType(42, 73));
    auto key2 = Traits::encode(PointType(42, 73));
    auto key3 = Traits::encode(PointType(73, 42));

    // Equality
    POLY_CHECK(key1 == key2);
    POLY_CHECK(!(key1 != key2));

    // Inequality
    POLY_CHECK(key1 != key3);
    POLY_CHECK(!(key1 == key3));

    cout << "Operator tests passed!" << endl;
  }

  cout << "\n=== All 2D Morton-key tests passed! ===" << endl;

  comm.finalize();
  return 0;
}
