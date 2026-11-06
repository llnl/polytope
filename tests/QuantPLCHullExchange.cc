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
    POLY_CHECK(record.second.second[0] == Q.quantize(Point2<double>(double(source), 0.0)));
  }

  comm.finalize();
  return 0;
}
