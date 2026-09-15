#include <filesystem>

namespace polytope {
//------------------------------------------------------------------------
// Generate testing variable data. Includes positions and ranks.
// TODO: Add materials
//------------------------------------------------------------------------
template<int Dimension, typename TessType>
void SiloWriterBase<Dimension, TessType>::generateTestVars() {
  const auto numCells = m_mesh.cells.size();
  std::vector<int> indices(numCells);
  std::array<std::vector<double>, Dimension> pos;
  for (auto& component : pos) component.resize(numCells);
  for (auto i = 0u; i < numCells; ++i) {
    for (int d = 0; d < Dimension; ++d) {
      pos[d][i] = m_mesh.points[i][d];
    }
    indices[i] = int(i);
  }
  std::array<std::string, 3> pos_names = {"x", "y", "z"};
  FieldMap<double> cellDoubleFields;
  for (int d = 0; d < Dimension; ++d) {
    cellDoubleFields[pos_names[d]] = pos[d];
  }
  FieldMap<int> cellIntFields;
  cellIntFields["cell_index"] = indices;
#ifdef POLYTOPE_ENABLE_MPI
  if (!m_mesh.cellRank.empty()) {
    std::vector<int> ranks = m_mesh.cellRank;
    cellIntFields["rank"] = ranks;
  } else {
    int rank = Communicator::getRank();
    std::vector<int> ranks(numCells, rank);
    cellIntFields["rank"] = ranks;
  }
#endif
  addField<double>(FieldCentering::Cell, cellDoubleFields);
  addField<int>(FieldCentering::Cell, cellIntFields);
}

//------------------------------------------------------------------------
// Determine the centering and write the appropriate field.
//------------------------------------------------------------------------
template<int Dimension, typename TessType>
template<typename FieldType>
void SiloWriterBase<Dimension, TessType>::
fieldWrite(const typename SiloWriterBase<Dimension, TessType>::FieldTypeMap<FieldType>& fields,
           const std::string& meshname,
           DBfile* file,
           DBoptlist* optlist) {
  const auto numFaces = m_mesh.faces.size();
  const auto numNodes = m_mesh.nodes.size();
  const auto numCells = m_mesh.cells.size();
  for (const auto& [centering, fieldmap] : fields) {
    const int centType = static_cast<int>(centering);
    if (centering == FieldCentering::Face ||
        centering == FieldCentering::Edge) {
      writeFields<FieldType>(fieldmap, meshname, file, numFaces, centType, optlist);
    } else if (centering == FieldCentering::Node) {
      writeFields<FieldType>(fieldmap, meshname, file, numNodes, centType, optlist);
    } else {
      writeFields<FieldType>(fieldmap, meshname, file, numCells, centType, optlist);
    }
  }
}

//------------------------------------------------------------------------
// Write the master file.
//------------------------------------------------------------------------
template<int Dimension, typename TessType>
void SiloWriterBase<Dimension, TessType>::
writeMasterFile(const std::string& masterFilename,
                const std::string& masterDirname,
                const std::vector<int>& ranksWithData,
                double& time,
                int cycle) {
  int coord_sys = DB_CARTESIAN;
  int driver = DB_HDF5;
  DBfile* file = nullptr;
  std::string ovlMFile = masterDirname + "/OvlTop.silo";
  if (m_overlinkType) {
    file = DBCreate(ovlMFile.c_str(), DB_CLOBBER, DB_LOCAL, "Master file", driver);
  } else {
    file = DBCreate(masterFilename.c_str(), DB_CLOBBER, DB_LOCAL, "Master file", driver);
  }

  // Build mesh names
  std::vector<std::string> procPaths = getProcPaths(masterDirname, ranksWithData);
  int nblocks = static_cast<int>(ranksWithData.size());
  const std::string meshname = getLocalMeshName();
  std::vector<char*> cellMeshNames;
  std::vector<char*> pointMeshNames;
  std::vector<int> varTypes(nblocks, DB_UCDVAR);
  std::vector<int> cellMeshTypes(nblocks, DB_UCDMESH);
  std::vector<int> pointMeshTypes(nblocks, DB_POINTMESH);
  for (const auto& p : procPaths) {
    cellMeshNames.push_back(strDup((p + meshname).c_str()));
    pointMeshNames.push_back(strDup((p + "points").c_str()));
  }
  DBoptlist* masteroptlist = DBMakeOptlist(10);
  if (cycle >= 0)
    DBAddOption(masteroptlist, DBOPT_CYCLE, &cycle);
  if (time >= 0.)
    DBAddOption(masteroptlist, DBOPT_DTIME, &time);
  std::string global_mesh_name = getGlobalMeshName();
  DBAddOption(masteroptlist, DBOPT_MMESH_NAME, global_mesh_name.data());
  DBAddOption(masteroptlist, DBOPT_COORDSYS, &coord_sys);

  DBPutMultimesh(file, global_mesh_name.c_str(), nblocks, cellMeshNames.data(), cellMeshTypes.data(), masteroptlist);
  DBPutMultimesh(file, "PPOINTS", nblocks, pointMeshNames.data(), pointMeshTypes.data(), masteroptlist);

  for (const auto& [centering, fieldmap] : m_doubleFields) {
    putCellVars(file, fieldmap, procPaths, nblocks, varTypes, masteroptlist);
  }
  for (const auto& [centering, fieldmap] : m_intFields) {
    putCellVars(file, fieldmap, procPaths, nblocks, varTypes, masteroptlist);
  }
#ifdef POLYTOPE_ENABLE_DEBUG
  std::vector<char*> nodeMeshNames;
  for (const auto& p : procPaths) {
    nodeMeshNames.push_back(strDup((p + "nodes").c_str()));
  }
  DBPutMultimesh(file, "NNODES", nblocks, nodeMeshNames.data(), pointMeshTypes.data(), masteroptlist);
  for (auto f = 0u; f < nodeMeshNames.size(); ++f) {
    free(nodeMeshNames[f]);
  }
#endif
  DBClose(file);
  // Clean up
  DBFreeOptlist(masteroptlist);
  for (int i = 0; i < nblocks; ++i) {
    free(cellMeshNames[i]);
    free(pointMeshNames[i]);
  }
  if (m_overlinkType) {
    std::error_code ec;
    std::filesystem::create_symlink(ovlMFile, masterFilename, ec);
  }
}
} // namespace polytope
