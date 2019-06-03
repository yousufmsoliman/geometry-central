#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/geometry_base.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <unordered_set>


using namespace geometrycentral;
using namespace geometrycentral::surface;
using std::cout;
using std::endl;

class HalfedgeGeometrySuite : public MeshAssetSuite {};

// ============================================================
// =============== Quantity management tests
// ============================================================

TEST_F(HalfedgeGeometrySuite, RefreshTest) {
  // TODO
}

TEST_F(HalfedgeGeometrySuite, PurgeTest) {
  // TODO
}

// ============================================================
// =============== Quantity tests
// ============================================================

// TODO change these to not explicitly construct a base object

TEST_F(HalfedgeGeometrySuite, VertexIndicesTest) {
  auto asset = getAsset("bob_small.ply");
  std::unique_ptr<HalfedgeMesh>& mesh = asset.mesh;
  cout << "building geom" << endl;
  GeometryBase geometry(*mesh);
  cout << "done building geom" << endl;


}
