#include "geometrycentral/surface/geometry_base.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/meshio.h"

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
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  GeometryBase& geometry = *asset.geometry;

  // Make sure the size is zero when empty
  EXPECT_EQ(geometry.vertexIndices.size(), 0);

  // Get them indices
  geometry.requireVertexIndices();
  EXPECT_EQ(geometry.vertexIndices.size(), mesh.nVertices());

  // Unrequire (but should not get rid of yet)
  geometry.unrequireVertexIndices();
  EXPECT_EQ(geometry.vertexIndices.size(), mesh.nVertices());


  // Purge actually deletes
  geometry.purgeQuantities();
  EXPECT_EQ(geometry.vertexIndices.size(), 0);
}

// ============================================================
// =============== Quantity tests
// ============================================================

// == Basic indices

TEST_F(HalfedgeGeometrySuite, VertexIndicesTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  GeometryBase& geometry = *asset.geometry;

  geometry.requireVertexIndices();
  for (Vertex v : mesh.vertices()) {
    EXPECT_GE(geometry.vertexIndices[v], 0);
    EXPECT_LT(geometry.vertexIndices[v], mesh.nVertices());
  }
}

TEST_F(HalfedgeGeometrySuite, HalfedgeIndicesTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  GeometryBase& geometry = *asset.geometry;

  geometry.requireHalfedgeIndices();
  for (Halfedge e : mesh.halfedges()) {
    EXPECT_GE(geometry.halfedgeIndices[e], 0);
    EXPECT_LT(geometry.halfedgeIndices[e], mesh.nHalfedges());
  }
}

TEST_F(HalfedgeGeometrySuite, CornerIndicesTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  GeometryBase& geometry = *asset.geometry;

  geometry.requireCornerIndices();
  for (Corner e : mesh.corners()) {
    EXPECT_GE(geometry.cornerIndices[e], 0);
    EXPECT_LT(geometry.cornerIndices[e], mesh.nCorners());
  }
}

TEST_F(HalfedgeGeometrySuite, EdgeIndicesTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  GeometryBase& geometry = *asset.geometry;

  geometry.requireEdgeIndices();
  for (Edge e : mesh.edges()) {
    EXPECT_GE(geometry.edgeIndices[e], 0);
    EXPECT_LT(geometry.edgeIndices[e], mesh.nEdges());
  }
}

TEST_F(HalfedgeGeometrySuite, FaceIndicesTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  GeometryBase& geometry = *asset.geometry;

  geometry.requireFaceIndices();
  for (Face e : mesh.faces()) {
    EXPECT_GE(geometry.faceIndices[e], 0);
    EXPECT_LT(geometry.faceIndices[e], mesh.nFaces());
  }
}

TEST_F(HalfedgeGeometrySuite, BoundaryLoopIndicesTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  GeometryBase& geometry = *asset.geometry;

  geometry.requireBoundaryLoopIndices();
  for (BoundaryLoop e : mesh.boundaryLoops()) {
    EXPECT_GE(geometry.boundaryLoopIndices[e], 0);
    EXPECT_LT(geometry.boundaryLoopIndices[e], mesh.nBoundaryLoops());
  }
}
