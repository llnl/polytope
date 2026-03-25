//------------------------------------------------------------------------------
// Unit tests of the polytope serialization methods.
//------------------------------------------------------------------------------
#include "polytope.hh"
#include "polytope_serialize.hh"
#include "polytope_test_utilities.hh"
#include "Point.hh"

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
checkSerialization(const T& val0) {
  vector<char> buffer;
  serialize(val0, buffer);
  T val1;
  vector<char>::const_iterator bufItr = buffer.begin();
  deserialize(val1, bufItr, buffer.end());
  // cout << "Check equivalence : " << val0 << " <=> " << val1 << endl;
  POLY_CHECK(bufItr == buffer.end());
  POLY_CHECK(val1 == val0);
}

//------------------------------------------------------------------------------
// The test itself.
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

  // Check out various primitive types we care about.
  checkSerialization((int) rand());
  checkSerialization((unsigned) abs(rand()));
  checkSerialization((uint32_t) abs(rand()));
  checkSerialization((uint64_t) abs(rand()));
  checkSerialization(random01());

  // Point2 types.
  checkSerialization(Point2<unsigned>(abs(rand()), abs(rand())));
  checkSerialization(Point2<uint32_t>((uint32_t) abs(rand()), (uint32_t) abs(rand())));
  checkSerialization(Point2<uint64_t>((uint64_t) abs(rand()), (uint64_t) abs(rand())));
  checkSerialization(Point2<double>(random01(), random01()));

  // Point3 types.
  checkSerialization(Point3<unsigned>(abs(rand()), abs(rand()), abs(rand())));
  checkSerialization(Point3<uint32_t>((uint32_t) abs(rand()), (uint32_t) abs(rand()), (uint32_t) abs(rand())));
  checkSerialization(Point3<uint64_t>((uint64_t) abs(rand()), (uint64_t) abs(rand()), (uint64_t) abs(rand())));
  checkSerialization(Point3<double>(random01(), random01(), random01()));

  //------------------------------------------------------------------------------
  // std::vector
  const size_t n1 = 100;
  {
    std::vector<int> val(n1);
    for (unsigned i = 0; i != n1; ++i) val[i] = rand();
    checkSerialization(val);
  }
  {
    std::vector<unsigned> val(n1);
    for (unsigned i = 0; i != n1; ++i) val[i] = abs(rand());
    checkSerialization(val);
  }
  {
    std::vector<uint32_t> val(n1);
    for (unsigned i = 0; i != n1; ++i) val[i] = abs(rand());
    checkSerialization(val);
  }
  {
    std::vector<uint64_t> val(n1);
    for (unsigned i = 0; i != n1; ++i) val[i] = abs(rand());
    checkSerialization(val);
  }
  {
    std::vector<double> val(n1);
    for (unsigned i = 0; i != n1; ++i) val[i] = numeric_limits<double>::max() * random01();
    checkSerialization(val);
  }
  {
    std::vector<Point2<uint64_t> > val(n1);
    for (unsigned i = 0; i != n1; ++i) val[i] = Point2<uint64_t>((uint64_t) abs(rand()), (uint64_t) abs(rand()));
    checkSerialization(val);
  }
  {
    std::vector<Point2<double> > val(n1);
    for (unsigned i = 0; i != n1; ++i) val[i] = Point2<double>(numeric_limits<double>::max() * random01(),
                                                               numeric_limits<double>::max() * random01());
    checkSerialization(val);
  }
  {
    std::vector<Point3<uint64_t> > val(n1);
    for (unsigned i = 0; i != n1; ++i) val[i] = Point3<uint64_t>((uint64_t) abs(rand()), (uint64_t) abs(rand()), (uint64_t) abs(rand()));
    checkSerialization(val);
  }
  {
    std::vector<Point3<double> > val(n1);
    for (unsigned i = 0; i != n1; ++i) val[i] = Point3<double>(numeric_limits<double>::max() * random01(),
                                                               numeric_limits<double>::max() * random01(),
                                                               numeric_limits<double>::max() * random01());
    checkSerialization(val);
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
    checkSerialization(val);
  }
  {
    std::vector<std::vector<unsigned> > val(n1, std::vector<unsigned>(n2));
    for (unsigned i = 0; i != n1; ++i) {
      for (unsigned j = 0; j != n2; ++j) {
        val[i][j] = abs(rand());
      }
    }
    checkSerialization(val);
  }
  {
    std::vector<std::vector<uint32_t> > val(n1, std::vector<uint32_t>(n2));
    for (unsigned i = 0; i != n1; ++i) {
      for (unsigned j = 0; j != n2; ++j) {
        val[i][j] = abs(rand());
      }
    }
    checkSerialization(val);
  }
  {
    std::vector<std::vector<uint64_t> > val(n1, std::vector<uint64_t>(n2));
    for (unsigned i = 0; i != n1; ++i) {
      for (unsigned j = 0; j != n2; ++j) {
        val[i][j] = abs(rand());
      }
    }
    checkSerialization(val);
  }
  {
    std::vector<std::vector<double> > val(n1, std::vector<double>(n2));
    for (unsigned i = 0; i != n1; ++i) {
      for (unsigned j = 0; j != n2; ++j) {
        val[i][j] = numeric_limits<double>::max() * random01();
      }
    }
    checkSerialization(val);
  }
  {
    std::vector<std::vector<Point2<uint64_t> > > val(n1, std::vector<Point2<uint64_t> >(n2));
    for (unsigned i = 0; i != n1; ++i) {
      for (unsigned j = 0; j != n2; ++j) {
        val[i][j] = Point2<uint64_t>((uint64_t) abs(rand()), (uint64_t) abs(rand()));
      }
    }
    checkSerialization(val);
  }
  {
    std::vector<std::vector<Point2<double> > > val(n1, std::vector<Point2<double> >(n2));
    for (unsigned i = 0; i != n1; ++i) {
      for (unsigned j = 0; j != n2; ++j) {
        val[i][j] = Point2<double>(numeric_limits<double>::max() * random01(),
                                   numeric_limits<double>::max() * random01());
      }
    }
    checkSerialization(val);
  }
  {
    std::vector<std::vector<Point3<uint64_t> > > val(n1, std::vector<Point3<uint64_t> >(n2));
    for (unsigned i = 0; i != n1; ++i) {
      for (unsigned j = 0; j != n2; ++j) {
        val[i][j] = Point3<uint64_t>((uint64_t) abs(rand()), (uint64_t) abs(rand()), (uint64_t) abs(rand()));
      }
    }
    checkSerialization(val);
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
    checkSerialization(val);
  }

  // That's all!
  cout << "PASS" << endl;
  return 0;
}
