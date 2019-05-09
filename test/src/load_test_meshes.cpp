#include "load_test_meshes.h"

#include "geometrycentral/mesh/polygon_soup_mesh.h"

using namespace geometrycentral;


std::unique_ptr<geometrycentral::HalfedgeMesh> mesh_from_soup(const std::vector<Vector3>& vertexPositions,
                                                              const std::vector<std::vector<size_t>>& faceIndices) {

  PolygonSoupMesh soup(faceIndices, vertexPositions);
  Geometry<Euclidean>* geom;
  std::unique_ptr<HalfedgeMesh> mesh(new HalfedgeMesh(soup, geom));
  delete geom;
  return mesh;
}



std::unique_ptr<HalfedgeMesh> load_tet() {
  // clang-format off
  std::vector<Vector3> vertexPositions = {
    Vector3{0., 0., 0.,},
    Vector3{1., 0., 0.,},
    Vector3{0., 1., 0.,},
    Vector3{0., 0., 1.,}
  };

  std::vector<std::vector<size_t>> faceIndices = {
    {0, 2, 1},
    {0, 1, 3},
    {0, 3, 2},
    {1, 2, 3}
  };
  // clang-format on

  return mesh_from_soup(vertexPositions, faceIndices);
}
