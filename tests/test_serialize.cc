//------------------------------------------------------------------------------
// Unit tests of the polytope serialization methods.
//------------------------------------------------------------------------------
#include "polytope.hh"
#include "polytope_serialize.hh"
#include "polytope_test_utilities.hh"
#include "Point.hh"
#include "HashKey.hh"

#include <vector>
#include <limits>

using namespace std;
using namespace polytope;

//------------------------------------------------------------------------------
// Templated helper method to reduce redundant code checking serializing
// various types.
//------------------------------------------------------------------------------
template<typename T>
void
checkSerialization(const T& val0, const string& description) {
  vector<char> buffer;
  serialize(val0, buffer);
  T val1;
  vector<char>::const_iterator bufItr = buffer.begin();
  deserialize(val1, bufItr, buffer.end());
  POLY_CHECK2(bufItr == buffer.end(),
              "Buffer not fully consumed after deserializing " << description);
  POLY_CHECK2(val1 == val0,
              "Deserialized value does not match original for " << description);
}

//------------------------------------------------------------------------------
// The test itself.
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

  // Check out various primitive types we care about.
  checkSerialization((int) rand(), "int");
  checkSerialization((unsigned) abs(rand()), "unsigned");
  checkSerialization((uint32_t) abs(rand()), "uint32_t");
  checkSerialization((uint64_t) abs(rand()), "uint64_t");
  checkSerialization(random01(), "double");

  // Point2 types.
  checkSerialization(Point2<unsigned>(abs(rand()), abs(rand())), "Point2<unsigned>");
  checkSerialization(Point2<uint32_t>((uint32_t) abs(rand()), (uint32_t) abs(rand())), "Point2<uint32_t>");
  checkSerialization(Point2<uint64_t>((uint64_t) abs(rand()), (uint64_t) abs(rand())), "Point2<uint64_t>");
  checkSerialization(Point2<double>(random01(), random01()), "Point2<double>");

  // Point3 types.
  checkSerialization(Point3<unsigned>(abs(rand()), abs(rand()), abs(rand())), "Point3<unsigned>");
  checkSerialization(Point3<uint32_t>((uint32_t) abs(rand()), (uint32_t) abs(rand()), (uint32_t) abs(rand())), "Point3<uint32_t>");
  checkSerialization(Point3<uint64_t>((uint64_t) abs(rand()), (uint64_t) abs(rand()), (uint64_t) abs(rand())), "Point3<uint64_t>");
  checkSerialization(Point3<double>(random01(), random01(), random01()), "Point3<double>");

  //------------------------------------------------------------------------------
  // HashKey2D types
  {
    // Test with small values
    HashKey2D key1;
    key1.interleave(Point2<uint32_t>(42, 73));
    checkSerialization(key1, "HashKey2D (small values)");

    // Test with large values
    HashKey2D key2;
    key2.interleave(Point2<uint32_t>(UINT32_MAX - 1, UINT32_MAX - 2));
    checkSerialization(key2, "HashKey2D (large values)");

    // Test with random values
    HashKey2D key3;
    key3.interleave(Point2<uint32_t>((uint32_t)abs(rand()), (uint32_t)abs(rand())));
    checkSerialization(key3, "HashKey2D (random values)");

    // Test with zero
    HashKey2D key4;
    key4.interleave(Point2<uint32_t>(0, 0));
    checkSerialization(key4, "HashKey2D (zero)");

    // Test with powers of 2
    HashKey2D key5;
    key5.interleave(Point2<uint32_t>(1024, 2048));
    checkSerialization(key5, "HashKey2D (powers of 2)");
  }

  //------------------------------------------------------------------------------
  // HashKey3D types
  {
    // Test with small values
    HashKey3D key1;
    key1.interleave(Point3<uint64_t>(42, 73, 101));
    checkSerialization(key1, "HashKey3D (small values)");

    // Test with large values (42-bit max)
    const uint64_t maxVal = (1ULL << 42) - 1;
    HashKey3D key2;
    key2.interleave(Point3<uint64_t>(maxVal, maxVal - 1, maxVal - 2));
    checkSerialization(key2, "HashKey3D (large 42-bit values)");

    // Test with random values
    HashKey3D key3;
    key3.interleave(Point3<uint64_t>((uint64_t)abs(rand()),
                                     (uint64_t)abs(rand()),
                                     (uint64_t)abs(rand())));
    checkSerialization(key3, "HashKey3D (random values)");

    // Test with zero
    HashKey3D key4;
    key4.interleave(Point3<uint64_t>(0, 0, 0));
    checkSerialization(key4, "HashKey3D (zero)");

    // Test with powers of 2
    HashKey3D key5;
    key5.interleave(Point3<uint64_t>(1024, 2048, 4096));
    checkSerialization(key5, "HashKey3D (powers of 2)");

    // Test with mixed high and low bits
    HashKey3D key6;
    key6.interleave(Point3<uint64_t>(1ULL << 40, 1ULL << 35, 1ULL << 30));
    checkSerialization(key6, "HashKey3D (mixed high/low bits)");
  }

  //------------------------------------------------------------------------------
  // std::vector
  const size_t n1 = 100;
  {
    std::vector<int> val(n1);
    for (unsigned i = 0; i != n1; ++i) val[i] = rand();
    checkSerialization(val, "std::vector<int>");
  }
  {
    std::vector<unsigned> val(n1);
    for (unsigned i = 0; i != n1; ++i) val[i] = abs(rand());
    checkSerialization(val, "std::vector<unsigned>");
  }
  {
    std::vector<uint32_t> val(n1);
    for (unsigned i = 0; i != n1; ++i) val[i] = abs(rand());
    checkSerialization(val, "std::vector<uint32_t>");
  }
  {
    std::vector<uint64_t> val(n1);
    for (unsigned i = 0; i != n1; ++i) val[i] = abs(rand());
    checkSerialization(val, "std::vector<uint64_t>");
  }
  {
    std::vector<double> val(n1);
    for (unsigned i = 0; i != n1; ++i) val[i] = numeric_limits<double>::max() * random01();
    checkSerialization(val, "std::vector<double>");
  }
  {
    std::vector<Point2<uint64_t> > val(n1);
    for (unsigned i = 0; i != n1; ++i) val[i] = Point2<uint64_t>((uint64_t) abs(rand()), (uint64_t) abs(rand()));
    checkSerialization(val, "std::vector<Point2<uint64_t>>");
  }
  {
    std::vector<Point2<double> > val(n1);
    for (unsigned i = 0; i != n1; ++i) val[i] = Point2<double>(numeric_limits<double>::max() * random01(),
                                                               numeric_limits<double>::max() * random01());
    checkSerialization(val, "std::vector<Point2<double>>");
  }
  {
    std::vector<Point3<uint64_t> > val(n1);
    for (unsigned i = 0; i != n1; ++i) val[i] = Point3<uint64_t>((uint64_t) abs(rand()), (uint64_t) abs(rand()), (uint64_t) abs(rand()));
    checkSerialization(val, "std::vector<Point3<uint64_t>>");
  }
  {
    std::vector<Point3<double> > val(n1);
    for (unsigned i = 0; i != n1; ++i) val[i] = Point3<double>(numeric_limits<double>::max() * random01(),
                                                               numeric_limits<double>::max() * random01(),
                                                               numeric_limits<double>::max() * random01());
    checkSerialization(val, "std::vector<Point3<double>>");
  }

  //------------------------------------------------------------------------------
  // std::vector<std::vector> >
  const size_t n2 = 10;
  {
    std::vector<std::vector<int> > val(n1, std::vector<int>(n2));
    for (unsigned i = 0; i != n1; ++i) {
      for (unsigned j = 0; j != n2; ++j) {
        val[i][j] = rand();
      }
    }
    checkSerialization(val, "std::vector<std::vector<int>>");
  }
  {
    std::vector<std::vector<unsigned> > val(n1, std::vector<unsigned>(n2));
    for (unsigned i = 0; i != n1; ++i) {
      for (unsigned j = 0; j != n2; ++j) {
        val[i][j] = abs(rand());
      }
    }
    checkSerialization(val, "std::vector<std::vector<unsigned>>");
  }
  {
    std::vector<std::vector<uint32_t> > val(n1, std::vector<uint32_t>(n2));
    for (unsigned i = 0; i != n1; ++i) {
      for (unsigned j = 0; j != n2; ++j) {
        val[i][j] = abs(rand());
      }
    }
    checkSerialization(val, "std::vector<std::vector<uint32_t>>");
  }
  {
    std::vector<std::vector<uint64_t> > val(n1, std::vector<uint64_t>(n2));
    for (unsigned i = 0; i != n1; ++i) {
      for (unsigned j = 0; j != n2; ++j) {
        val[i][j] = abs(rand());
      }
    }
    checkSerialization(val, "std::vector<std::vector<uint64_t>>");
  }
  {
    std::vector<std::vector<double> > val(n1, std::vector<double>(n2));
    for (unsigned i = 0; i != n1; ++i) {
      for (unsigned j = 0; j != n2; ++j) {
        val[i][j] = numeric_limits<double>::max() * random01();
      }
    }
    checkSerialization(val, "std::vector<std::vector<double>>");
  }
  {
    std::vector<std::vector<Point2<uint64_t> > > val(n1, std::vector<Point2<uint64_t> >(n2));
    for (unsigned i = 0; i != n1; ++i) {
      for (unsigned j = 0; j != n2; ++j) {
        val[i][j] = Point2<uint64_t>((uint64_t) abs(rand()), (uint64_t) abs(rand()));
      }
    }
    checkSerialization(val, "std::vector<std::vector<Point2<uint64_t>>>");
  }
  {
    std::vector<std::vector<Point2<double> > > val(n1, std::vector<Point2<double> >(n2));
    for (unsigned i = 0; i != n1; ++i) {
      for (unsigned j = 0; j != n2; ++j) {
        val[i][j] = Point2<double>(numeric_limits<double>::max() * random01(),
                                   numeric_limits<double>::max() * random01());
      }
    }
    checkSerialization(val, "std::vector<std::vector<Point2<double>>>");
  }
  {
    std::vector<std::vector<Point3<uint64_t> > > val(n1, std::vector<Point3<uint64_t> >(n2));
    for (unsigned i = 0; i != n1; ++i) {
      for (unsigned j = 0; j != n2; ++j) {
        val[i][j] = Point3<uint64_t>((uint64_t) abs(rand()), (uint64_t) abs(rand()), (uint64_t) abs(rand()));
      }
    }
    checkSerialization(val, "std::vector<std::vector<Point3<uint64_t>>>");
  }
  {
    std::vector<std::vector<Point3<double> > > val(n1, std::vector<Point3<double> >(n2));
    for (unsigned i = 0; i != n1; ++i) {
      for (unsigned j = 0; j != n2; ++j) {
        val[i][j] = Point3<double>(numeric_limits<double>::max() * random01(),
                                   numeric_limits<double>::max() * random01(),
                                   numeric_limits<double>::max() * random01());
      }
    }
    checkSerialization(val, "std::vector<std::vector<Point3<double>>>");
  }

  //------------------------------------------------------------------------------
  // std::vector<HashKey2D>
  {
    std::vector<HashKey2D> val(n1);
    for (unsigned i = 0; i != n1; ++i) {
      HashKey2D key;
      key.interleave(Point2<uint32_t>((uint32_t)abs(rand()), (uint32_t)abs(rand())));
      val[i] = key;
    }
    checkSerialization(val, "std::vector<HashKey2D>");
  }

  //------------------------------------------------------------------------------
  // std::vector<HashKey3D>
  {
    std::vector<HashKey3D> val(n1);
    for (unsigned i = 0; i != n1; ++i) {
      HashKey3D key;
      key.interleave(Point3<uint64_t>((uint64_t)abs(rand()),
                                      (uint64_t)abs(rand()),
                                      (uint64_t)abs(rand())));
      val[i] = key;
    }
    checkSerialization(val, "std::vector<HashKey3D>");
  }

  //------------------------------------------------------------------------------
  // std::vector<std::vector<HashKey2D>>
  {
    std::vector<std::vector<HashKey2D>> val(n2);
    for (unsigned i = 0; i != n2; ++i) {
      val[i].resize(n2);
      for (unsigned j = 0; j != n2; ++j) {
        HashKey2D key;
        key.interleave(Point2<uint32_t>((uint32_t)abs(rand()), (uint32_t)abs(rand())));
        val[i][j] = key;
      }
    }
    checkSerialization(val, "std::vector<std::vector<HashKey2D>>");
  }

  //------------------------------------------------------------------------------
  // std::vector<std::vector<HashKey3D>>
  {
    std::vector<std::vector<HashKey3D>> val(n2);
    for (unsigned i = 0; i != n2; ++i) {
      val[i].resize(n2);
      for (unsigned j = 0; j != n2; ++j) {
        HashKey3D key;
        key.interleave(Point3<uint64_t>((uint64_t)abs(rand()),
                                        (uint64_t)abs(rand()),
                                        (uint64_t)abs(rand())));
        val[i][j] = key;
      }
    }
    checkSerialization(val, "std::vector<std::vector<HashKey3D>>");
  }

  // That's all!
  
  return 0;
}
