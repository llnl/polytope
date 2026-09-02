// test_Circumcenter

#include "polytope.hh"
#include "polytope_test_utilities.hh"
#include "Communicator.hh"
#include "GeomUtils.hh"
#include <vector>

using namespace polytope;

namespace {
void createTriangle(const double loval,
                    const double hival,
                    std::array<Point2<double>, 3>& tri1,
                    std::array<Point2<double>, 3>& tri2) {
  Point2<double> ad(loval, loval);
  Point2<double> bd(hival, loval);
  Point2<double> cd(loval, hival);
  Point2<double> dd(hival, hival);
  tri1[0] = ad;
  tri1[1] = bd;
  tri1[2] = cd;
  tri2[0] = dd;
  tri2[1] = cd;
  tri2[2] = bd;
}

template<typename IntType>
Point2<IntType> circumcenter2d(std::array<Point2<double>, 3> tri) {
  auto a = tri[0];
  auto b = tri[1];
  auto c = tri[2];
  return round<2, IntType>(circumcenter(a, b, c));
}

template<typename IntType>
Point3<IntType> circumcenter3d(std::array<Point3<double>, 4> tri) {
  auto a = tri[0];
  auto b = tri[1];
  auto c = tri[2];
  auto d = tri[3];
  return round<3, IntType>(circumcenter(a, b, c, d));
}

void createPrism(const double loval,
                 const double hival,
                 const double height,
                 std::array<Point3<double>, 4>& tri1,
                 std::array<Point3<double>, 4>& tri2) {
  Point3<double> hix(hival, loval, loval);
  Point3<double> hiy(loval, hival, loval);
  Point3<double> hiz(loval, loval, hival);
  // We must compute an origin to ensure the centroid if the tetrahedron
  // is at the centroid of the triangle formed from hix, hiy, and hiz
  Point3<double> centroid = triangleCentroid(hix, hiy, hiz);
  Point3<double> norm = normal(hix, hiy, hiz);
  std::cout << "norm " << norm << std::endl;
  Point3<long double> p = static_cast<long double>(4.)*centroid.template type_cast<long double>()
    - hix.template type_cast<long double>()
    - hiy.template type_cast<long double>()
    - hiz.template type_cast<long double>();
  Point3<double> origin = p.template type_cast<double>() + height*norm;
  Point3<double> origin2 = origin - 2.*height*norm;
  tri1[0] = origin;
  tri1[1] = hix;
  tri1[2] = hiy;
  tri1[3] = hiz;
  tri2[0] = hix;
  tri2[1] = hiy;
  tri2[2] = hiz;
  tri2[3] = origin2;
}
};



// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv) {
  using IntType2D = QuantizedCoordinate<2>;
  using IntType3D = QuantizedCoordinate<3>;
  auto& comm = Communicator::instance();
  comm.init(argc, argv);

#ifdef POLYTOPE_ENABLE_TRIANGLE
  {
    cout << "\nTest 2D circumcenter\n";
    Point2<double> min(-0.5, -0.5);
    Point2<double> max(0.5, 0.5);
    auto& q = Quantizer<2>::instance();
    q.init(min, max, 0.1);
    double loval = -0.375;
    double hival = -0.25;
    std::array<Point2<double>, 3> tri1, tri2;
    createTriangle(loval, hival, tri1, tri2);
    std::array<Point2<double>, 3> qtri1, qtri2;
    for (int i = 0; i < 3; ++i) {
      qtri1[i] = q.quantize(tri1[i]).template type_cast<double>();
      qtri2[i] = q.quantize(tri2[i]).template type_cast<double>();
    }
    auto c1 = circumcenter2d<IntType2D>(qtri1);
    auto c2 = circumcenter2d<IntType2D>(qtri2);
    POLY_CHECK2(c1 == c2, "Circumcenters differ between the triangles " << c1 << " " << c2);
  }
#endif

#ifdef POLYTOPE_ENABLE_TETGEN
  {
    cout << "\nTest 3D circumcenter\n";
    Point3<double> min(-0.5, -0.5, -0.5);
    Point3<double> max(0.5, 0.5, 0.5);
    auto& q = Quantizer<3>::instance();
    q.init(min, max, 0.1);
    double loval = 0.;//-0.375;
    double hival = 0.25;//-0.25;
    double height = 0.9;
    std::array<Point3<double>, 4> tri1, tri2;
    createPrism(loval, hival, height, tri1, tri2);
    std::array<Point3<double>, 4> qtri1, qtri2;
    for (int i = 0; i < 4; ++i) {
      std::cout << tri1[i] << std::endl;
      std::cout << tri2[i] << std::endl;
      qtri1[i] = q.quantize(tri1[i]).template type_cast<double>();
      qtri2[i] = q.quantize(tri2[i]).template type_cast<double>();
    }
    auto c1 = circumcenter3d<IntType3D>(qtri1);
    auto c2 = circumcenter3d<IntType3D>(qtri2);
    std::cout << c1 << " " << c2 << std::endl;
  }
#endif

  comm.finalize();
  return 0;
}
