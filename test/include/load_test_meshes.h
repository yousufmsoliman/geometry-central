#pragma once

#include "geometrycentral/surface/geometry.h"
#include "geometrycentral/surface/halfedge_mesh.h"

#include "gtest/gtest.h"

// A mesh used for testing
struct MeshAsset {
  MeshAsset(std::string localPath);

  std::string name = "Unnamed_Mesh_Asset";
  std::string sourcePath = "unknown";
  std::unique_ptr<geometrycentral::surface::HalfedgeMesh> mesh;
  std::unique_ptr<geometrycentral::surface::Geometry<geometrycentral::surface::Euclidean>> geometry;
  bool hasBoundary = false;
  bool isTriangular = true;
  bool isSimplicialComplex = true;

  MeshAsset copy();
};

// Loads test meshes from disk
class MeshAssetSuite : public ::testing::Test {
protected:
  static void SetUpTestSuite();
  // static void TearDownTestSuite();
  // virtual void SetUp();
  // virtual void TearDown();

  static std::vector<MeshAsset> allMeshAssets;


	// Get various groups for meshes
	
};
