#include "geometrycentral/mesh/halfedge_mesh.h"

namespace geometrycentral {

HalfedgeMeshDataTransfer::HalfedgeMeshDataTransfer(HalfedgeMesh* oldMesh_, HalfedgeMesh* newMesh_)
    : oldMesh(oldMesh_), newMesh(newMesh_) {
  heMap = HalfedgeData<Halfedge>(oldMesh);
  vMap = VertexData<Vertex>(oldMesh);
  eMap = EdgeData<Edge>(oldMesh);
  fMap = FaceData<Face>(oldMesh);
}

HalfedgeMeshDataTransfer::HalfedgeMeshDataTransfer() : oldMesh(nullptr), newMesh(nullptr) {}


void HalfedgeMeshDataTransfer::generateReverseMaps() {

  // Halfedges
  heMapBack = HalfedgeData<Halfedge>(newMesh);
  for (Halfedge he : oldMesh->allHalfedges()) {
    heMapBack[heMap[he]] = he;
  }

  // Vertices
  vMapBack = VertexData<Vertex>(newMesh);
  for (Vertex v : oldMesh->vertices()) {
    vMapBack[vMap[v]] = v;
  }

  // Edges
  eMapBack = EdgeData<Edge>(newMesh);
  for (Edge e : oldMesh->edges()) {
    eMapBack[eMap[e]] = e;
  }

  // Faces
  fMapBack = FaceData<Face>(newMesh);
  for (Face f : oldMesh->faces()) {
    fMapBack[fMap[f]] = f;
  }
  for (Face f : oldMesh->boundaryLoops()) {
    fMapBack[fMap[f]] = f;
  }
}

} // namespace geometrycentral
