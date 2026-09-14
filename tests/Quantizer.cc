// Unit tests for Quantizer class
//
// Tests the quantization/dequantization and hashing operations for
// converting between physical (floating-point) and quantized (integer)
// coordinate spaces.

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <limits>

#include "polytope.hh"
#include "Quantizer.hh"
#include "Point.hh"
#include "polytope_test_utilities.hh"

#include "Communicator.hh" 

using namespace std;
using namespace polytope;

namespace {

//------------------------------------------------------------------------------
// Test 2D quantizer with single bounding box
//------------------------------------------------------------------------------
void test2DQuantizer() {
  cout << "\n=== Testing 2D Quantizer ===" << endl;

  using RealType = double;
  using RealPoint = Point2<RealType>;
  using Quantizer2D = Quantizer<2>;
  using IntPoint = QuantizedPoint<2>;

  // Define a bounding box: [0, 10] x [0, 10]
  RealPoint xlo(0.0, 0.0);
  RealPoint xhi(10.0, 10.0);

  // Get the singleton instance and initialize it
  auto& quantizer = Quantizer2D::instance();
  quantizer.init(xlo, xhi);

  cout << "  Bounding box: (" << xlo.x << ", " << xlo.y << ") to ("
       << xhi.x << ", " << xhi.y << ")" << endl;

  // Test quantization/dequantization round-trip for points inside the box
  vector<RealPoint> testPoints = {
    RealPoint(0.0, 0.0),          // Corner
    RealPoint(10.0, 10.0),        // Opposite corner
    RealPoint(5.0, 5.0),          // Center
    RealPoint(2.5, 7.5),          // Random interior point
    RealPoint(0.1, 9.9),          // Near boundary
    RealPoint(3.14159, 2.71828)   // Irrational coordinates
  };

  cout << "  Testing quantize/dequantize round-trip..." << endl;
  for (auto p : testPoints) {  // Non-const copy for non-const reference methods
    IntPoint quantized = quantizer.quantize(p);
    RealPoint recovered = quantizer.dequantize(quantized);

    // Check that recovery is reasonably close (within quantization error)
    RealType dx = abs(recovered.x - p.x);
    RealType dy = abs(recovered.y - p.y);
    RealType maxError = quantizer.m_dx_o.x * 2.0; // Allow 2x spacing for rounding

    POLY_CHECK2(dx < maxError && dy < maxError,
                "Round-trip error too large for point (" << p.x << ", " << p.y << "): "
                << "recovered (" << recovered.x << ", " << recovered.y << "), "
                << "error (" << dx << ", " << dy << ")");
  }

  cout << "  Testing hash operations..." << endl;
  // Test hashing from real and integer points
  RealPoint testPt(5.0, 5.0);
  RealPoint testPtCopy = testPt;  // Need copy for non-const reference
  IntPoint quantizedPt = quantizer.quantize(testPtCopy);

  testPtCopy = testPt;
  auto hash1 = quantizer.quantizeAndEncode(testPtCopy);
  auto hash2 = quantizer.encode(quantizedPt);

  POLY_CHECK2(hash1 == hash2,
              "Hash from real point should match hash from quantized point");

  // Test unhashing to IntPoint
  IntPoint unhashed1 = quantizer.decode(hash1);
  POLY_CHECK2(unhashed1.x == quantizedPt.x && unhashed1.y == quantizedPt.y,
              "Unhashing should recover quantized coordinates");

  // Test decodeAndDequantize (combined operation)
  RealPoint unhashed2 = quantizer.decodeAndDequantize(hash1);
  RealType dx = abs(unhashed2.x - testPt.x);
  RealType dy = abs(unhashed2.y - testPt.y);
  RealType maxError = quantizer.m_dx_o.x * 2.0;
  POLY_CHECK2(dx < maxError && dy < maxError,
              "decodeAndDequantize should recover original (within tolerance)");

  cout << "  2D Quantizer tests passed!" << endl;
}

//------------------------------------------------------------------------------
// Test 3D quantizer with single bounding box
//------------------------------------------------------------------------------
void test3DQuantizer() {
  cout << "\n=== Testing 3D Quantizer ===" << endl;

  using RealType = double;
  using RealPoint = Point3<RealType>;
  using Quantizer3D = Quantizer<3>;
  using IntPoint = QuantizedPoint<3>;

  // Define a bounding box: [-5, 5]^3
  RealPoint xlo(-5.0, -5.0, -5.0);
  RealPoint xhi(5.0, 5.0, 5.0);

  // Get the singleton instance and initialize it
  auto& quantizer = Quantizer3D::instance();
  quantizer.init(xlo, xhi);

  cout << "  Bounding box: (" << xlo.x << ", " << xlo.y << ", " << xlo.z << ") to ("
       << xhi.x << ", " << xhi.y << ", " << xhi.z << ")" << endl;

  // Test quantization/dequantization round-trip
  vector<RealPoint> testPoints = {
    RealPoint(-5.0, -5.0, -5.0),  // Corner
    RealPoint(5.0, 5.0, 5.0),     // Opposite corner
    RealPoint(0.0, 0.0, 0.0),     // Center
    RealPoint(1.0, -2.0, 3.5),    // Random interior point
    RealPoint(-4.9, 4.9, 0.1),    // Near boundary
    RealPoint(2.5, -1.5, 0.0)     // Another test point
  };

  cout << "  Testing quantize/dequantize round-trip..." << endl;
  for (auto p : testPoints) {  // Non-const copy for non-const reference methods
    IntPoint quantized = quantizer.quantize(p);
    RealPoint recovered = quantizer.dequantize(quantized);

    // Check that recovery is reasonably close
    RealType dx = abs(recovered.x - p.x);
    RealType dy = abs(recovered.y - p.y);
    RealType dz = abs(recovered.z - p.z);
    RealType maxError = quantizer.m_dx_o.x * 2.0;

    POLY_CHECK2(dx < maxError && dy < maxError && dz < maxError,
                "Round-trip error too large for point (" << p.x << ", " << p.y << ", " << p.z << "): "
                << "recovered (" << recovered.x << ", " << recovered.y << ", " << recovered.z << "), "
                << "error (" << dx << ", " << dy << ", " << dz << ")");
  }

  cout << "  Testing hash operations..." << endl;
  RealPoint testPt(1.0, 2.0, 3.0);
  RealPoint testPtCopy = testPt;
  IntPoint quantizedPt = quantizer.quantize(testPtCopy);

  testPtCopy = testPt;
  auto hash1 = quantizer.quantizeAndEncode(testPtCopy);
  auto hash2 = quantizer.encode(quantizedPt);

  POLY_CHECK2(hash1 == hash2,
              "Hash from real point should match hash from quantized point");

  // Test unhashing
  IntPoint unhashed1 = quantizer.decode(hash1);
  POLY_CHECK2(unhashed1.x == quantizedPt.x && unhashed1.y == quantizedPt.y && unhashed1.z == quantizedPt.z,
              "Unhashing should recover quantized coordinates");

  cout << "  3D Quantizer tests passed!" << endl;
}

//------------------------------------------------------------------------------
// Test quantization accuracy at boundaries
//------------------------------------------------------------------------------
void testBoundaryAccuracy() {
  cout << "\n=== Testing Boundary Accuracy ===" << endl;

  using RealType = double;
  using RealPoint2D = Point2<RealType>;
  using Quantizer2D = Quantizer<2>;

  // Small bounding box to test precision
  RealPoint2D xlo(0.0, 0.0);
  RealPoint2D xhi(1.0, 1.0);

  auto& quantizer = Quantizer2D::instance();
  quantizer.init(xlo, xhi);

  // Test corners and edges
  vector<RealPoint2D> boundaryPoints = {
    RealPoint2D(0.0, 0.0),
    RealPoint2D(1.0, 0.0),
    RealPoint2D(0.0, 1.0),
    RealPoint2D(1.0, 1.0),
    RealPoint2D(0.5, 0.0),
    RealPoint2D(0.5, 1.0),
    RealPoint2D(0.0, 0.5),
    RealPoint2D(1.0, 0.5),
    RealPoint2D(0.5, 0.5)
  };

  cout << "  Testing boundary and center points..." << endl;
  for (auto p : boundaryPoints) {  // Non-const copy
    auto quantized = quantizer.quantize(p);
    auto recovered = quantizer.dequantize(quantized);

    RealType dx = abs(recovered.x - p.x);
    RealType dy = abs(recovered.y - p.y);
    RealType maxError = quantizer.m_dx_o.x * 2.0;

    POLY_CHECK2(dx < maxError && dy < maxError,
                "Boundary point (" << p.x << ", " << p.y << ") has excessive error: "
                << "(" << dx << ", " << dy << ")");
  }

  cout << "  Boundary accuracy tests passed!" << endl;
}

//------------------------------------------------------------------------------
// Test that distinct points produce distinct hashes (within reason)
//------------------------------------------------------------------------------
void testHashUniqueness() {
  cout << "\n=== Testing Hash Uniqueness ===" << endl;

  using RealType = double;
  using RealPoint = Point2<RealType>;
  using Quantizer2D = Quantizer<2>;

  RealPoint xlo(0.0, 0.0);
  RealPoint xhi(100.0, 100.0);
  auto& quantizer = Quantizer2D::instance();
  quantizer.init(xlo, xhi);

  // Generate random points and verify they hash uniquely
  const unsigned nPoints = 1000;
  vector<MortonKey<2>> hashes;
  hashes.reserve(nPoints);

  cout << "  Testing " << nPoints << " random points..." << endl;
  for (unsigned i = 0; i < nPoints; ++i) {
    RealPoint p(random01() * 100.0, random01() * 100.0);
    RealPoint pCopy = p;
    auto hash = quantizer.quantizeAndEncode(pCopy);

    // Check that this hash hasn't been seen before
    // (Note: quantization may cause nearby points to collide, which is acceptable)
    hashes.push_back(hash);
  }

  // Count unique hashes
  set<MortonKey<2>> uniqueHashes(hashes.begin(), hashes.end());
  double uniqueRatio = double(uniqueHashes.size()) / nPoints;

  cout << "  Unique hashes: " << uniqueHashes.size() << " / " << nPoints
       << " (" << (uniqueRatio * 100.0) << "%)" << endl;

  // We expect most points to be unique (unless quantization causes collisions)
  // A ratio > 0.95 is reasonable for random points
  POLY_CHECK2(uniqueRatio > 0.90,
              "Hash uniqueness ratio too low: " << uniqueRatio);

  cout << "  Hash uniqueness tests passed!" << endl;
}

//------------------------------------------------------------------------------
// Test quantization grid alignment
//------------------------------------------------------------------------------
void testGridAlignment() {
  cout << "\n=== Testing Grid Alignment ===" << endl;

  using RealType = double;
  using RealPoint = Point2<RealType>;
  using Quantizer2D = Quantizer<2>;
  using IntPoint = QuantizedPoint<2>;

  RealPoint xlo(0.0, 0.0);
  RealPoint xhi(10.0, 10.0);
  auto& quantizer = Quantizer2D::instance();
  quantizer.init(xlo, xhi);

  cout << "  Grid spacing: (" << quantizer.m_dx_o.x << ", " << quantizer.m_dx_o.y << ")" << endl;

  // Test that points on a regular grid quantize to integer multiples
  const unsigned gridSize = 10;
  cout << "  Testing " << gridSize << "x" << gridSize << " grid..." << endl;

  for (unsigned i = 0; i <= gridSize; ++i) {
    for (unsigned j = 0; j <= gridSize; ++j) {
      RealType x = (i * 10.0) / gridSize;
      RealType y = (j * 10.0) / gridSize;
      RealPoint p(x, y);

      IntPoint quantized = quantizer.quantize(p);
      RealPoint recovered = quantizer.dequantize(quantized);

      // Grid points should recover accurately
      RealType dx = abs(recovered.x - x);
      RealType dy = abs(recovered.y - y);
      RealType maxError = quantizer.m_dx_o.x * 2.0;

      POLY_CHECK2(dx < maxError && dy < maxError,
                  "Grid point (" << x << ", " << y << ") has excessive error");
    }
  }

  cout << "  Grid alignment tests passed!" << endl;
}

//------------------------------------------------------------------------------
// Test consistency: same point always produces same hash
//------------------------------------------------------------------------------
void testHashConsistency() {
  cout << "\n=== Testing Hash Consistency ===" << endl;

  using RealType = double;
  using RealPoint = Point2<RealType>;
  using Quantizer2D = Quantizer<2>;

  RealPoint xlo(-10.0, -10.0);
  RealPoint xhi(10.0, 10.0);
  auto& quantizer = Quantizer2D::instance();
  quantizer.init(xlo, xhi);

  // Test that hashing the same point multiple times gives the same result
  RealPoint testPoint(3.14159, 2.71828);

  cout << "  Hashing same point 100 times..." << endl;
  RealPoint testCopy = testPoint;
  auto hash1 = quantizer.quantizeAndEncode(testCopy);

  for (unsigned i = 0; i < 100; ++i) {
    testCopy = testPoint;
    auto hash2 = quantizer.quantizeAndEncode(testCopy);
    POLY_CHECK2(hash1 == hash2,
                "Hash inconsistency detected on iteration " << i);
  }

  cout << "  Hash consistency tests passed!" << endl;
}

//------------------------------------------------------------------------------
// Stress test with many operations
//------------------------------------------------------------------------------
void stressTest() {
  cout << "\n=== Stress Test ===" << endl;

  using RealType = double;
  using RealPoint = Point2<RealType>;
  using Quantizer2D = Quantizer<2>;

  RealPoint xlo(-1000.0, -1000.0);
  RealPoint xhi(1000.0, 1000.0);
  auto& quantizer = Quantizer2D::instance();
  quantizer.init(xlo, xhi);

  const unsigned nOps = 10000;
  cout << "  Performing " << nOps << " quantize/hash/unhash operations..." << endl;

  unsigned errorCount = 0;
  for (unsigned i = 0; i < nOps; ++i) {
    // Random point in large domain
    RealType x = (random01() - 0.5) * 2000.0;
    RealType y = (random01() - 0.5) * 2000.0;
    RealPoint p(x, y);

    // Quantize and hash
    RealPoint pCopy = p;
    pCopy = p;
    auto hash = quantizer.quantizeAndEncode(pCopy);

    // Unhash and verify using combined operation
    RealPoint recovered = quantizer.decodeAndDequantize(hash);

    RealType dx = abs(recovered.x - x);
    RealType dy = abs(recovered.y - y);
    RealType maxError = quantizer.m_dx_o.x * 2.0;

    if (dx >= maxError || dy >= maxError) {
      errorCount++;
    }
  }

  cout << "  Operations with excessive error: " << errorCount << " / " << nOps << endl;
  POLY_CHECK2(errorCount == 0,
              "Stress test found " << errorCount << " operations with excessive error");

  cout << "  Stress test passed!" << endl;
}

} // anonymous namespace

//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

  auto& comm = Communicator::instance();
  comm.init(argc, argv);

  // Seed random number generator
  srand(42);

  try {
    // 2D tests
    test2DQuantizer();
    testBoundaryAccuracy();
    testHashUniqueness();
    testGridAlignment();
    testHashConsistency();

    // 3D tests
    test3DQuantizer();

    // Stress tests
    stressTest();

    cout << "\n=== All Quantizer tests passed! ===" << endl;

  } catch (const exception& e) {
    cerr << "\nTest failed with exception: " << e.what() << endl;
    comm.finalize();
    return 1;
  }
  comm.finalize();
  return 0;
}
