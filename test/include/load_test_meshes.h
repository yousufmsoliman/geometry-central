#pragma once

#include "geometrycentral/surface/halfedge_mesh.h"


std::unique_ptr<geometrycentral::surface::HalfedgeMesh>
mesh_from_soup(const std::vector<geometrycentral::Vector3>& vertexPositions,
               const std::vector<std::vector<size_t>>& faceIndices);

std::unique_ptr<geometrycentral::surface::HalfedgeMesh> load_tet();
