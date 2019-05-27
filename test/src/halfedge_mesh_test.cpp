
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/meshio.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <iostream>
#include <string>


using namespace geometrycentral;
using namespace geometrycentral::surface;
using std::cout;
using std::endl;

class HalfedgeMeshSuite : public MeshAssetSuite {};

TEST_F(HalfedgeMeshSuite, LoadMeshTest) {
  // std::unique_ptr<HalfedgeMesh> m = load_tet();
  // m->validateConnectivity();
}


TEST(HalfedgeBasics, ValidateClosedMeshTest) {
  // std::unique_ptr<HalfedgeMesh> m = load_tet();
  // m->validateConnectivity();
}
