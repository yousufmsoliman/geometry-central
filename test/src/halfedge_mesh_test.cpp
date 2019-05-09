#include <iostream>
#include <string>

#include "gtest/gtest.h"

#include "load_test_meshes.h"

#include "geometrycentral/mesh/halfedge_mesh.h"

using namespace geometrycentral;
using std::cout;
using std::endl;

TEST(HalfedgeBasics, ValidateClosedMeshTest) {
  std::unique_ptr<HalfedgeMesh> m = load_tet();
  m->validateConnectivity();
}
