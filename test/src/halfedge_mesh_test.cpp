
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

class HalfedgeMeshSuite : public MeshAssetSuite {};

// ============================================================
// =============== Basic validation tests
// ============================================================

TEST_F(HalfedgeMeshSuite, ValidateClosedMeshTest) {
  for (MeshAsset& a : closedMeshes()) {
    a.printThyName();
    a.mesh->validateConnectivity();
  }
}

TEST_F(HalfedgeMeshSuite, ValidateBoundaryMeshTest) {
  for (MeshAsset& a : boundaryMeshes()) {
    a.printThyName();
    a.mesh->validateConnectivity();
  }
}

// ============================================================
// =============== Basic iterator tests
// ============================================================

TEST_F(HalfedgeMeshSuite, IterateVerticesTest) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    std::unordered_set<Vertex> seenElements;
    for (Vertex e : a.mesh->vertices()) {
      seenElements.insert(e);
    }

    EXPECT_EQ(seenElements.size(), a.mesh->nVertices());
  }
}

TEST_F(HalfedgeMeshSuite, IterateHalfedgesTest) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    std::unordered_set<Halfedge> seenElements;
    for (Halfedge e : a.mesh->halfedges()) {
      seenElements.insert(e);
    }

    EXPECT_EQ(seenElements.size(), a.mesh->nHalfedges());
  }
}

TEST_F(HalfedgeMeshSuite, IterateInteriorHalfedgesTest) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    std::unordered_set<Halfedge> seenElements;
    for (Halfedge e : a.mesh->interiorHalfedges()) {
      seenElements.insert(e);
    }

    EXPECT_EQ(seenElements.size(), a.mesh->nInteriorHalfedges());
  }
}

TEST_F(HalfedgeMeshSuite, IterateExteriorHalfedgesTest) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    std::unordered_set<Halfedge> seenElements;
    for (Halfedge e : a.mesh->exteriorHalfedges()) {
      seenElements.insert(e);
    }

    EXPECT_EQ(seenElements.size(), a.mesh->nExteriorHalfedges());
  }
}

TEST_F(HalfedgeMeshSuite, IterateCornersTest) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    std::unordered_set<Corner> seenElements;
    for (Corner e : a.mesh->corners()) {
      seenElements.insert(e);
    }

    EXPECT_EQ(seenElements.size(), a.mesh->nCorners());
  }
}

TEST_F(HalfedgeMeshSuite, IterateEdgesTest) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    std::unordered_set<Edge> seenElements;
    for (Edge e : a.mesh->edges()) {
      seenElements.insert(e);
    }

    EXPECT_EQ(seenElements.size(), a.mesh->nEdges());
  }
}

TEST_F(HalfedgeMeshSuite, IterateFacesTest) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    std::unordered_set<Face> seenElements;
    for (Face e : a.mesh->faces()) {
      seenElements.insert(e);
    }

    EXPECT_EQ(seenElements.size(), a.mesh->nFaces());
  }
}

TEST_F(HalfedgeMeshSuite, IterateBoundaryLoopsTest) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    std::unordered_set<BoundaryLoop> seenElements;
    for (BoundaryLoop e : a.mesh->boundaryLoops()) {
      seenElements.insert(e);
    }

    EXPECT_EQ(seenElements.size(), a.mesh->nBoundaryLoops());
  }
}
