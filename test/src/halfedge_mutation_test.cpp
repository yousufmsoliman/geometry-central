
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

// Split a few edges on a bunch of meshes
TEST_F(HalfedgeMutationSuite, InsertVertexAlongEdgeTest) {

  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    int count = 10;
    int indInc = static_cast<int>(std::ceil(a.mesh->nVertices() / static_cast<double>(count)));

    int ind = 0;
    for (int i = 0; i < count; i++) {

      // Insert along an edge
      Edge e = a.mesh->edge(ind);
      a.mesh->insertVertexAlongEdge(e);
      a.mesh->validateConnectivity();

      ind = (ind + indInc) % a.mesh->nVertices();
    }
  }
}


// Insert a vertex along every edge and triangulate (not-quite subdivision)
TEST_F(HalfedgeMutationSuite, InsertVertexAndTriangulateSubdivideTest) {

  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    // Split every edge
    std::vector<Edge> origEdges;
    for (Edge e : a.mesh->edges()) {
      origEdges.push_back(e);
    }
    for (Edge e : origEdges) {
      a.mesh->insertVertexAlongEdge(e);
    }

    a.mesh->validateConnectivity();

    // Triangulate
    // TODO this loops while modifying. Do we allow that?
    for (Face f : a.mesh->faces()) {
      a.mesh->triangulate(f);
    }

    a.mesh->validateConnectivity();
  }
}

// Split every edge and then flip (regular subdivision)
TEST_F(HalfedgeMutationSuite, SplitFlipSubdivie) {

  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    VertexData<char> isNewVertex(*a.mesh, false);
    for (Vertex v : a.mesh->vertices()) {
      isNewVertex[v] = true;
    }

    // Split every edge
    std::vector<Edge> origEdges;
    for (Edge e : a.mesh->edges()) {
      origEdges.push_back(e);
    }
    for (Edge e : origEdges) {
      a.mesh->splitEdge(e);
    }
    a.mesh->validateConnectivity();

    // Flip edges between old and new
    for (Edge e : a.mesh->edges()) {
      if (isNewVertex[e.halfedge().vertex()] != isNewVertex[e.halfedge().twin().vertex()]) {
        a.mesh->flip(e);
      }
    }
    a.mesh->validateConnectivity();

    // Should yield subdivision
    for (Face f : a.mesh->faces()) {
      EXPECT_TRUE(f.isTriangle());
    }
  }
}

// Split a few edges on a bunch of meshes
TEST_F(HalfedgeMutationSuite, EdgeSplitTest) {

  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    int count = 10;
    int indInc = static_cast<int>(std::ceil(a.mesh->nVertices() / static_cast<double>(count)));

    int splitInd = 0;
    for (int i = 0; i < count; i++) {

      // Split an edge
      Edge eSplit = a.mesh->edge(splitInd);
      a.mesh->splitEdge(eSplit);
      a.mesh->validateConnectivity();

      splitInd = (splitInd + indInc) % a.mesh->nVertices();
    }
  }
}
