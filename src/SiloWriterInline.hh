#include <filesystem>

namespace polytope {
//------------------------------------------------------------------------
// Generate testing variable data. Includes positions and ranks.
//------------------------------------------------------------------------
template<int Dimension, typename TessType>
void SiloWriterBase<Dimension, TessType>::generateTestVars() {
  const auto numCells = m_mesh.cells.size();
  std::vector<int> indices(numCells);
  std::array<std::vector<double>, Dimension> pos;
  for (auto& component : pos) component.resize(numCells);
  for (auto i = 0u; i < numCells; ++i) {
    for (int d = 0; d < Dimension; ++d) {
      pos[d][i] = static_cast<double>(m_mesh.points[i][d]);
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
  std::vector<std::string> matNames(1, "O2");
  std::vector<std::vector<double>> matVF(numCells, std::vector<double>(1, 1.));
  addMaterials(matNames, matVF);
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
// Configure the non-default Overlink attributes for a field.
//------------------------------------------------------------------------
template<int Dimension, typename TessType>
void SiloWriterBase<Dimension, TessType>::
setOverlinkAttributes(const std::string& fieldName,
                      int scaling,
                      int linking,
                      int metric) {
  if ((scaling != OverlinkAttr.at("ATTR_INTENSIVE") &&
       scaling != OverlinkAttr.at("ATTR_EXTENSIVE")) ||
      (linking != OverlinkAttr.at("ATTR_FIRST_ORDER") &&
       linking != OverlinkAttr.at("ATTR_SECOND_ORDER")) ||
      metric < 0 || metric > 2) {
    throw std::invalid_argument("Invalid Overlink field attributes");
  }
  m_overlinkFieldAttributes[fieldName] = {scaling, linking, metric};
}

//------------------------------------------------------------------------
// Write Overlink's per-variable remapping metadata to the master file.
//------------------------------------------------------------------------
template<int Dimension, typename TessType>
void SiloWriterBase<Dimension, TessType>::
writeOverlinkAttributes(DBfile* file, DBoptlist* optlist) {
  std::vector<const char*> variableNames;
  std::vector<int> attributeLengths;
  std::vector<int> attributes;

  const auto appendFields = [this, &variableNames, &attributeLengths, &attributes]
    (const auto& fields, int dataType) {
      for (const auto& [centering, fieldmap] : fields) {
        // Overlink VAR_ATTRIBUTES supports nodal and zonal variables only.
        if (centering != FieldCentering::Node &&
            centering != FieldCentering::Cell) continue;
        const int overlinkCentering =
          centering == FieldCentering::Node ?
          OverlinkAttr.at("ATTR_NODAL") : OverlinkAttr.at("ATTR_ZONAL");
        for (const auto& field : fieldmap) {
          const auto& name = field.first;
          const auto attrIter = m_overlinkFieldAttributes.find(name);
          const auto scaling = attrIter == m_overlinkFieldAttributes.end() ?
            OverlinkAttr.at("ATTR_INTENSIVE") : attrIter->second.scaling;
          const auto linking = attrIter == m_overlinkFieldAttributes.end() ?
            OverlinkAttr.at("ATTR_SECOND_ORDER") : attrIter->second.linking;
          const auto metric = attrIter == m_overlinkFieldAttributes.end() ?
            0 : attrIter->second.metric;
          variableNames.push_back(name.c_str());
          attributeLengths.push_back(6);
          attributes.insert(attributes.end(),
                            {overlinkCentering, scaling, linking, 0, dataType, metric});
        }
      }
    };

  appendFields(m_doubleFields, OverlinkAttr.at("ATTR_FLOAT"));
  appendFields(m_intFields, OverlinkAttr.at("ATTR_INTEGER"));
  if (!variableNames.empty()) {
    DBPutCompoundarray(file,
                       "VAR_ATTRIBUTES",
                       variableNames.data(),
                       attributeLengths.data(),
                       static_cast<int>(variableNames.size()),
                       attributes.data(),
                       static_cast<int>(attributes.size()),
                       DB_INT,
                       optlist);
  }
}

//------------------------------------------------------------------------
// Convert and save the material volume fractions into Silo appropriate forms
// matVF is each cells volume fractions stored as matVF[cell][material]
//------------------------------------------------------------------------
template<int Dimension, typename TessType>
void SiloWriterBase<Dimension, TessType>::
addMaterials(const std::vector<std::string>& matNames,
             const std::vector<std::vector<double>>& matVF) {
  const auto nmat = matNames.size();
  const auto nzones = m_mesh.cells.size();
  m_matNames = matNames;
  POLY_CHECK(nzones == matVF.size());
  m_matlist.clear();
  m_mix_next.clear();
  m_mix_mat.clear();
  m_mix_zone.clear();
  m_mix_vf.clear();
  m_matlist.resize(nzones);
  // Loop over each cell
  for (auto i = 0u; i < nzones; ++i) {
    const auto& fractions = matVF[i];
    POLY_CHECK(fractions.size() == nmat);
    double sum = 0.;
    std::vector<std::pair<size_t, double>> activeMaterials;
    // Loop over each material
    for (auto j = 0u; j < nmat; ++j) {
      const double fraction = fractions[j];
      POLY_ASSERT(fraction >= 0.);
      if (fraction > m_matTolerance) {
        activeMaterials.emplace_back(j, fraction);
        sum += fraction;
      }
    }
    POLY_ASSERT(!activeMaterials.empty());
    POLY_ASSERT(std::abs(sum - 1.) < m_matTolerance);
    // Check if it is a clean (non-mixed) zone
    if (activeMaterials.size() == 1) {
      m_matlist[i] = activeMaterials.front().first + 1;
      continue;
    }
    // Mixed zone. Silo uses a negative value to identify the first
    // mixed-material entry. The index is 1-based in the negative value.
    m_matlist[i] = -static_cast<int>(m_mix_mat.size()) - 1;
    int previousMixEntry = -1;
    for (const auto& [material, fraction] : activeMaterials) {
      const int mixEntry =
        static_cast<int>(m_mix_mat.size());
      m_mix_mat.push_back(material + 1);
      m_mix_vf.push_back(fraction);
      m_mix_zone.push_back(i);
      m_mix_next.push_back(0);
      if (previousMixEntry >= 0) {
        m_mix_next[previousMixEntry] = mixEntry + 1;
      }
      previousMixEntry = mixEntry;
    }
  }
}

//------------------------------------------------------------------------
// Write the material data.
//------------------------------------------------------------------------
template<int Dimension, typename TessType>
void SiloWriterBase<Dimension, TessType>::
writeMaterialsToFile(const std::string& meshname,
                     const std::string& matname,
                     DBfile* file,
                     DBoptlist* optlist) {
  const auto nzones = m_mesh.cells.size();
  const auto nmat = m_matNames.size();
  if (nmat == 0) {
    return;
  }
  const auto mixlen = m_mix_vf.size();
  std::vector<int> matnos(nmat);
  for (int i = 0; i < nmat; ++i) {
    matnos[i] = i + 1;
  }
  std::vector<int> dims(1, nzones);
  std::vector<char*> matNamesPtr;
  for (auto& name : m_matNames) {
    matNamesPtr.push_back(name.data());
  }
  DBAddOption(optlist, DBOPT_MATNAMES, matNamesPtr.data());
  DBPutMaterial(file,
                matname.c_str(),
                meshname.c_str(),
                nmat,
                matnos.data(),
                m_matlist.data(),
                dims.data(),
                1,
                mixlen == 0 ? nullptr : m_mix_next.data(),
                mixlen == 0 ? nullptr : m_mix_mat.data(),
                mixlen == 0 ? nullptr : m_mix_zone.data(),
                mixlen == 0 ? nullptr : m_mix_vf.data(),
                mixlen,
                DB_DOUBLE,
                optlist);
  DBClearOption(optlist, DBOPT_MATNAMES);
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
  const std::string matname = getLocalMatName();
  std::vector<char*> cellMeshNames;
  std::vector<char*> pointMeshNames;
  std::vector<char*> materialMeshNames;
  std::vector<int> varTypes(nblocks, DB_UCDVAR);
  std::vector<int> cellMeshTypes(nblocks, DB_UCDMESH);
  std::vector<int> pointMeshTypes(nblocks, DB_POINTMESH);
  for (const auto& p : procPaths) {
    cellMeshNames.push_back(strDup((p + meshname).c_str()));
    pointMeshNames.push_back(strDup((p + "points").c_str()));
    materialMeshNames.push_back(strDup((p + matname).c_str()));
  }
  DBoptlist* masteroptlist = DBMakeOptlist(10);
  if (cycle >= 0)
    DBAddOption(masteroptlist, DBOPT_CYCLE, &cycle);
  if (time >= 0.)
    DBAddOption(masteroptlist, DBOPT_DTIME, &time);
  std::string global_mesh_name = getGlobalMeshName();
  std::string global_mat_name = getGlobalMatName();
  DBAddOption(masteroptlist, DBOPT_MMESH_NAME, global_mesh_name.data());
  DBAddOption(masteroptlist, DBOPT_COORDSYS, &coord_sys);

  DBPutMultimesh(file, global_mesh_name.c_str(), nblocks, cellMeshNames.data(), cellMeshTypes.data(), masteroptlist);
  DBPutMultimesh(file, "PPOINTS", nblocks, pointMeshNames.data(), pointMeshTypes.data(), masteroptlist);

  // Add field data
  for (const auto& [centering, fieldmap] : m_doubleFields) {
    putCellVars(file, fieldmap, procPaths, nblocks, varTypes, masteroptlist);
  }
  for (const auto& [centering, fieldmap] : m_intFields) {
    putCellVars(file, fieldmap, procPaths, nblocks, varTypes, masteroptlist);
  }
  if (m_overlinkType) {
    writeOverlinkAttributes(file, masteroptlist);
  }
  // Add material data
  const auto nmat = m_matNames.size();
  if (nmat > 0) {
    std::vector<int> matnos(nmat);
    std::vector<char*> matNamesPtr;
    for (auto i = 0u; i < nmat; ++i) {
      matnos[i] = i + 1;
      auto& name = m_matNames[i];
      matNamesPtr.push_back(name.data());
    }
    int nmatnos = static_cast<int>(nmat);
    DBAddOption(masteroptlist, DBOPT_NMATNOS, &nmatnos);
    DBAddOption(masteroptlist, DBOPT_MATNOS, matnos.data());
    DBAddOption(masteroptlist, DBOPT_MATNAMES, matNamesPtr.data());
    DBPutMultimat(file, global_mat_name.c_str(), nblocks, materialMeshNames.data(), masteroptlist);
    DBClearOption(masteroptlist, DBOPT_MATNAMES);
    DBClearOption(masteroptlist, DBOPT_MATNOS);
    DBClearOption(masteroptlist, DBOPT_NMATNOS);
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
    free(materialMeshNames[i]);
  }
  if (m_overlinkType) {
    std::error_code ec;
    std::filesystem::create_symlink(ovlMFile, masterFilename, ec);
  }
}
} // namespace polytope
