
#include "Communicator.hh"
#include "QuantTessellation.hh"
#include "Tessellator.hh"

#ifdef POLYTOPE_ENABLE_BOOST
#include "BoostTessellator.hh"
#endif
#ifdef POLYTOPE_ENABLE_TRIANGLE
#include "TriangleTessellator.hh"
#endif

#include <iostream>
#include <memory>
#include <string>

using namespace polytope;

namespace {

std::unique_ptr<Tessellator<2, double>>
makeTessellator(const std::string& name) {
#ifdef POLYTOPE_ENABLE_BOOST
  if (name == "BoostTessellator") {
    return std::unique_ptr<Tessellator<2, double>>(new BoostTessellator());
  }
#endif

#ifdef POLYTOPE_ENABLE_TRIANGLE
  if (name == "TriangleTessellator") {
    return std::unique_ptr<Tessellator<2, double>>(new TriangleTessellator());
  }
#endif

  return std::unique_ptr<Tessellator<2, double>>();
}

}

int main(int argc, char** argv) {
  auto& comm = Communicator::instance();
  comm.init(argc, argv);

  if (argc < 2) {
    std::cout << "Must provide an escape pod file" << std::endl;
    Communicator::abort();
  }
  std::string inputFile = argv[1];
  auto& q = Quantizer<2>::instance();
  q.init(Point2<double>(0.), Point2<double>(1.));

  QuantTessellation<2> quantMesh;
  QuantPLC<2> qplc;
  std::string tessellatorName;
  quantMesh.loadEscapePod(inputFile, qplc, tessellatorName);
  quantMesh.m_isEscapePod = true;

  if (tessellatorName.empty()) {
    std::cout << "Escape pod does not specify a serial tessellator" << std::endl;
    Communicator::abort();
  }

  auto tessellator = makeTessellator(tessellatorName);
  if (!tessellator) {
    std::cout << "Unsupported or unavailable serial tessellator '"
              << tessellatorName << "'" << std::endl;
    Communicator::abort();
  }

  tessellator->tessellateQuantized(quantMesh);
  if (!qplc.empty()) {
    quantMesh.clipTessellation(qplc, *tessellator);
  } else {
    std::cout << "Escape pod does not contain QPLC data; skipping clipping" << std::endl;
  }
#ifdef POLYTOPE_ENABLE_SILO
  std::string prefix = "escapefile";
  std::map<std::string, std::vector<double>> cellFields;
  size_t meshSize = quantMesh.cells.size();
  std::vector<double> index(meshSize);
  std::vector<double> genx (meshSize);
  std::vector<double> geny (meshSize);
  for (int i = 0; i < meshSize; ++i) {
    index[i] = double(i);
    genx[i] = quantMesh.points[i].x;
    geny[i] = quantMesh.points[i].y;
  }
  cellFields["cell_index"] = index;
  cellFields["gen_x"     ] = genx;
  cellFields["gen_y"     ] = geny;
  std::map<int, std::map<std::string, std::vector<double>>> fields;
  fields[DB_ZONECENT] = cellFields;
  SiloWriter<2, QuantTessellation<2>>::write(quantMesh, fields, prefix, 1, 0., 1);
#endif
  comm.finalize();
  return 0;
}
