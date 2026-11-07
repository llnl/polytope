#include <vector>

#include "Communicator.hh"
#include "ParallelUtils.hh"
#include "polytope_test_utilities.hh"

using namespace polytope;

int main(int argc, char** argv) {
  auto& comm = Communicator::instance();
  comm.init(argc, argv);

  const int rank = Communicator::getRank();
  const int nranks = Communicator::getNRanks();
  Quantizer<2> Q;
  Q.init(Point2<double>(-1.0, -1.0), Point2<double>(nranks + 2.0, 2.0));
  const std::vector<Point2<double>> realPoints = {
    {double(rank) + 0.1, 0.1},
    {double(rank), 0.0}, {double(rank) + 0.5, 0.0},
    {double(rank), 0.5}
  };
  auto points = QuantPLC<2>::quantizePoints(Q, realPoints);
  QuantPLC<2> hull(Q);
  hull.makeConvex(points);

  const auto gathered = allGatherHulls(true, hull, points);
  POLY_CHECK(int(gathered.size()) == nranks);
  for (int source = 0; source < nranks; ++source) {
    const auto& record = gathered[source];
    POLY_CHECK(record.first);
    POLY_CHECK(record.second.second.size() == 3);
    POLY_CHECK(record.second.first.isValid(record.second.second));
    POLY_CHECK(record.second.first.m_convex);
    const std::vector<Point2<double>> sourceRealPoints = {
      {double(source) + 0.1, 0.1},
      {double(source), 0.0}, {double(source) + 0.5, 0.0},
      {double(source), 0.5}
    };
    const auto sourcePoints = QuantPLC<2>::quantizePoints(Q, sourceRealPoints);
    QuantPLC<2> sourceHull(Q);
    sourceHull.makeConvex(sourcePoints);
    POLY_CHECK(record.second.first.getCell(record.second.second) ==
               sourceHull.getCell(sourcePoints));
  }

  Quantizer<3> Q3;
  Q3.init(Point3<double>(-1.0, -1.0, -1.0),
          Point3<double>(nranks + 2.0, 2.0, 2.0));
  const std::vector<Point3<double>> realPoints3 = {
    {double(rank) + 0.1, 0.1, 0.1},
    {double(rank), 0.0, 0.0}, {double(rank) + 0.5, 0.0, 0.0},
    {double(rank), 0.5, 0.0}, {double(rank), 0.0, 0.5}
  };
  auto points3 = QuantPLC<3>::quantizePoints(Q3, realPoints3);
  QuantPLC<3> hull3(Q3);
  hull3.facets = {{1, 3, 2}, {1, 2, 4}, {2, 3, 4}, {3, 1, 4}};
  hull3.m_convex = true;

  const auto gathered3 = allGatherHulls(true, hull3, points3);
  POLY_CHECK(int(gathered3.size()) == nranks);
  for (int source = 0; source < nranks; ++source) {
    const auto& record = gathered3[source];
    POLY_CHECK(record.first);
    POLY_CHECK(record.second.second.size() == 4);
    const std::vector<Point3<double>> sourceRealPoints = {
      {double(source) + 0.1, 0.1, 0.1},
      {double(source), 0.0, 0.0}, {double(source) + 0.5, 0.0, 0.0},
      {double(source), 0.5, 0.0}, {double(source), 0.0, 0.5}
    };
    const auto sourcePoints = QuantPLC<3>::quantizePoints(Q3, sourceRealPoints);
    POLY_CHECK(record.second.first.getCell(record.second.second) ==
               hull3.getCell(sourcePoints));
  }

  comm.finalize();
  return 0;
}
