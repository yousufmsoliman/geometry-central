
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/meshio.h"

#include "geometrycentral/surface/base_geometry_interface.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/surface/extrinsic_geometry_interface.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <unordered_set>


using namespace geometrycentral;
using namespace geometrycentral::surface;
using std::cout;
using std::endl;

class HalfedgeMutationSuite : public MeshAssetSuite {};

// Flip a few edges on a bunch of meshes
TEST_F(HalfedgeMutationSuite, EdgeFlipTest) {

  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    int count = 10;
    int indInc = static_cast<int>(std::ceil(a.mesh->nVertices() / static_cast<double>(count)));

    int flipInd = 0;
    for (int i = 0; i < count; i++) {

      // Flip an edge
      Edge eFlip = a.mesh->edge(flipInd);
      a.mesh->flip(eFlip);
      a.mesh->validateConnectivity();

      flipInd = (flipInd + indInc) % a.mesh->nVertices();
    }
  }
}

// Flip a lot of edges on one mesh without boundary
TEST_F(HalfedgeMutationSuite, EdgeFlipClosedManyTest) {

  auto asset = getAsset("sphere_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;

  int count = 1000;
  int indInc = static_cast<int>(std::ceil(mesh.nVertices() / static_cast<double>(count)));

  int flipInd = 0;
  for (int i = 0; i < count; i++) {

    // Flip an edge
    Edge eFlip = mesh.edge(flipInd);
    mesh.flip(eFlip);
    mesh.validateConnectivity();

    flipInd = (flipInd + 1) % mesh.nVertices();
  }
}
