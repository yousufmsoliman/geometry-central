#include "load_test_meshes.h"

#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/polygon_soup_mesh.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

using std::cout;
using std::endl;

// Helpers
namespace {
std::string guessNiceNameFromPath(std::string fullname) {
  size_t startInd = 0;
  for (std::string sep : {"/", "\\"}) {
    size_t pos = fullname.rfind(sep);
    if (pos != std::string::npos) {
      startInd = std::max(startInd, pos + 1);
    }
  }

  return fullname.substr(startInd, fullname.size() - startInd);
};
} // namespace

MeshAsset::MeshAsset(std::string localPath) {
  name = guessNiceNameFromPath(localPath);
  std::string fullPath = std::string(GC_TEST_ASSETS_ABS_PATH) + "/" + name;
  cout << "  -- info: Loading mesh asset " << name << " from " << fullPath << endl;
  sourcePath = fullPath;
  std::tie(mesh, geometry) = loadMesh(fullPath);
}

// Static storage for mesh assets
std::vector<MeshAsset> MeshAssetSuite::allMeshAssets;

void MeshAssetSuite::SetUpTestSuite() {
  allMeshAssets.emplace_back("tet.obj");
  allMeshAssets.emplace_back("spot.ply");
  allMeshAssets.emplace_back("bob_small.ply");
  allMeshAssets.emplace_back("lego.ply");
}

/*
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

  return std::unique_ptr<HalfedgeMesh>(new HalfedgeMesh(faceIndices));
}
*/
