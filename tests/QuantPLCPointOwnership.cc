#include <vector>

#include "Communicator.hh"
#include "QuantPLC.hh"
#include "polytope_test_utilities.hh"

using namespace polytope;

int main(int argc, char** argv) {
  auto& comm = Communicator::instance();
  comm.init(argc, argv);

  Quantizer<2> Q;
  Q.init(Point2<double>(-1.0, -1.0), Point2<double>(5.0, 5.0));
  const std::vector<Point2<double>> realPoints = {
    {0.0, 0.0}, {4.0, 0.0}, {4.0, 4.0}, {0.0, 4.0},
    {1.0, 1.0}, {1.0, 3.0}, {3.0, 3.0}, {3.0, 1.0},
    {4.5, 4.5}
  };
  auto points = QuantPLC<2>::quantizePoints(Q, realPoints);
  PLC<2> plc;
  plc.facets = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
  plc.holes = {{{4, 5}, {5, 6}, {6, 7}, {7, 4}}};
  QuantPLC<2> qplc(Q, points, plc);

  POLY_CHECK(points.size() == 9);
  POLY_CHECK(qplc.within(Point2<double>(0.5, 0.5), points));
  POLY_CHECK(!qplc.within(Point2<double>(2.0, 2.0), points));
  POLY_CHECK(qplc.getHolePoints(points).size() == 1);

  qplc.reduce(points);
  POLY_CHECK(points.size() == 8);
  POLY_CHECK(qplc.m_reduced);
  POLY_CHECK(qplc.within(Point2<double>(0.5, 0.5), points));
  POLY_CHECK(!qplc.within(Point2<double>(2.0, 2.0), points));

  auto hull = QuantPLC<2>(Q);
  hull.makeConvex(points);
  POLY_CHECK(hull.isValid(points));
  POLY_CHECK(hull.m_convex);
  POLY_CHECK(QuantPLC<2>::convexPLCIntersection(hull, points, hull, points));

  comm.finalize();
  return 0;
}
