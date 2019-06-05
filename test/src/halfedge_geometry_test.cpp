
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
  BaseGeometryInterface& geometry = *asset.geometry;

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

// Simple tests which ensure that the quantity can be computed and is a reasonable range.

// == Basic indices

TEST_F(HalfedgeGeometrySuite, VertexIndicesTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  BaseGeometryInterface& geometry = *asset.geometry;

  geometry.requireVertexIndices();
  for (Vertex v : mesh.vertices()) {
    EXPECT_GE(geometry.vertexIndices[v], 0);
    EXPECT_LT(geometry.vertexIndices[v], mesh.nVertices());
  }
}

TEST_F(HalfedgeGeometrySuite, HalfedgeIndicesTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  BaseGeometryInterface& geometry = *asset.geometry;

  geometry.requireHalfedgeIndices();
  for (Halfedge e : mesh.halfedges()) {
    EXPECT_GE(geometry.halfedgeIndices[e], 0);
    EXPECT_LT(geometry.halfedgeIndices[e], mesh.nHalfedges());
  }
}

TEST_F(HalfedgeGeometrySuite, CornerIndicesTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  BaseGeometryInterface& geometry = *asset.geometry;

  geometry.requireCornerIndices();
  for (Corner e : mesh.corners()) {
    EXPECT_GE(geometry.cornerIndices[e], 0);
    EXPECT_LT(geometry.cornerIndices[e], mesh.nCorners());
  }
}

TEST_F(HalfedgeGeometrySuite, EdgeIndicesTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  BaseGeometryInterface& geometry = *asset.geometry;

  geometry.requireEdgeIndices();
  for (Edge e : mesh.edges()) {
    EXPECT_GE(geometry.edgeIndices[e], 0);
    EXPECT_LT(geometry.edgeIndices[e], mesh.nEdges());
  }
}

TEST_F(HalfedgeGeometrySuite, FaceIndicesTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  BaseGeometryInterface& geometry = *asset.geometry;

  geometry.requireFaceIndices();
  for (Face e : mesh.faces()) {
    EXPECT_GE(geometry.faceIndices[e], 0);
    EXPECT_LT(geometry.faceIndices[e], mesh.nFaces());
  }
}

TEST_F(HalfedgeGeometrySuite, BoundaryLoopIndicesTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  BaseGeometryInterface& geometry = *asset.geometry;

  geometry.requireBoundaryLoopIndices();
  for (BoundaryLoop e : mesh.boundaryLoops()) {
    EXPECT_GE(geometry.boundaryLoopIndices[e], 0);
    EXPECT_LT(geometry.boundaryLoopIndices[e], mesh.nBoundaryLoops());
  }
}

// == Intrinsic geometry

TEST_F(HalfedgeGeometrySuite, EdgeLengthsTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  IntrinsicGeometryInterface& geometry = *asset.geometry;

  geometry.requireEdgeLengths();
  for (Edge e : mesh.edges()) {
    EXPECT_GT(geometry.edgeLengths[e], 0);
    EXPECT_TRUE(std::isfinite(geometry.edgeLengths[e]));
  }
}


TEST_F(HalfedgeGeometrySuite, FaceAreasTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  IntrinsicGeometryInterface& geometry = *asset.geometry;

  geometry.requireFaceAreas();
  for (Face e : mesh.faces()) {
    EXPECT_GE(geometry.faceAreas[e], 0);
    EXPECT_TRUE(std::isfinite(geometry.faceAreas[e]));
  }
}

TEST_F(HalfedgeGeometrySuite, VertexDualAreasTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  IntrinsicGeometryInterface& geometry = *asset.geometry;

  geometry.requireVertexDualAreas();
  for (Vertex e : mesh.vertices()) {
    EXPECT_GE(geometry.vertexDualAreas[e], 0);
    EXPECT_TRUE(std::isfinite(geometry.vertexDualAreas[e]));
  }
}


TEST_F(HalfedgeGeometrySuite, CornerAnglesTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  IntrinsicGeometryInterface& geometry = *asset.geometry;

  geometry.requireCornerAngles();
  for (Corner e : mesh.corners()) {
    EXPECT_GE(geometry.cornerAngles[e], 0);
    EXPECT_TRUE(std::isfinite(geometry.cornerAngles[e]));
  }
}

TEST_F(HalfedgeGeometrySuite, VertexAngleSumsTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  IntrinsicGeometryInterface& geometry = *asset.geometry;

  geometry.requireVertexAngleSums();
  for (Vertex e : mesh.vertices()) {
    EXPECT_GE(geometry.vertexAngleSums[e], 0);
    EXPECT_TRUE(std::isfinite(geometry.vertexAngleSums[e]));
  }
}


TEST_F(HalfedgeGeometrySuite, CornerScaledAnglesTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  IntrinsicGeometryInterface& geometry = *asset.geometry;

  geometry.requireCornerScaledAngles();
  for (Corner e : mesh.corners()) {
    EXPECT_GE(geometry.cornerScaledAngles[e], 0);
    EXPECT_TRUE(std::isfinite(geometry.cornerScaledAngles[e]));
  }
}


TEST_F(HalfedgeGeometrySuite, VertexGaussianCurvaturesTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  IntrinsicGeometryInterface& geometry = *asset.geometry;

  geometry.requireVertexGaussianCurvatures();
  for (Vertex e : mesh.vertices()) {
    EXPECT_TRUE(std::isfinite(geometry.vertexGaussianCurvatures[e]));
  }
}


TEST_F(HalfedgeGeometrySuite, FaceGaussianCurvaturesTest) {
  auto asset = getAsset("bob_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  IntrinsicGeometryInterface& geometry = *asset.geometry;

  geometry.requireFaceGaussianCurvatures();
  for (Face e : mesh.faces()) {
    EXPECT_TRUE(std::isfinite(geometry.faceGaussianCurvatures[e]));
  }
}

TEST_F(HalfedgeGeometrySuite, HalfedgeCotanWeightsTest) {
  auto asset = getAsset("lego.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  IntrinsicGeometryInterface& geometry = *asset.geometry;

  geometry.requireHalfedgeCotanWeights();
  for (Halfedge e : mesh.halfedges()) {
    EXPECT_TRUE(std::isfinite(geometry.halfedgeCotanWeights[e]));
  }
}


TEST_F(HalfedgeGeometrySuite, EdgeCotanWeightsTest) {
  auto asset = getAsset("lego.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  IntrinsicGeometryInterface& geometry = *asset.geometry;

  geometry.requireEdgeCotanWeights();
  for (Edge e : mesh.edges()) {
    EXPECT_TRUE(std::isfinite(geometry.edgeCotanWeights[e]));
  }
}


// ============================================================
// =============== Geometry tests
// ============================================================

// More interesting geometry tests which check invariants, etc


// Check that the vertex curvatures return the value expected by Gauss-Bonnet
TEST_F(HalfedgeGeometrySuite, VertexGaussianCurvaturesSumTest) {
  for (auto& asset : closedMeshes()) {
    if (!asset.isTriangular) continue;

    asset.printThyName();
    HalfedgeMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireVertexGaussianCurvatures();
    double curvatureSum = 0.;
    for (Vertex e : mesh.vertices()) {
      curvatureSum += geometry.vertexGaussianCurvatures[e];
    }

    double gaussBonnetCurvature = 2. * PI * mesh.eulerCharacteristic();

    EXPECT_NEAR(curvatureSum, gaussBonnetCurvature, 1e-2);
  }
}

// Check that the face curvatures return the value expected by Gauss-Bonnet
TEST_F(HalfedgeGeometrySuite, FaceGaussianCurvaturesSumTest) {
  for (auto& asset : closedMeshes()) {
    if (!asset.isTriangular) continue;

    asset.printThyName();
    HalfedgeMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireFaceGaussianCurvatures();
    double curvatureSum = 0.;
    for (Face f : mesh.faces()) {
      curvatureSum += geometry.faceGaussianCurvatures[f];
    }

    double gaussBonnetCurvature = 2. * PI * mesh.eulerCharacteristic();

    EXPECT_NEAR(curvatureSum, gaussBonnetCurvature, 1e-2);
  }
}
