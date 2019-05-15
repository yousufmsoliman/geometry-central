#include "geometrycentral/mesh/halfedge_mesh.h"

namespace geometrycentral {

template <class T>
VertexData<T> HalfedgeMeshDataTransfer::transfer(VertexData<T>& inData) {

  VertexData<T> outData(newMesh);
  for (Vertex v : oldMesh->vertices()) {
    outData[vMap[v]] = inData[v];
  }

  return outData;
}

template <class T>
FaceData<T> HalfedgeMeshDataTransfer::transfer(FaceData<T>& inData) {

  FaceData<T> outData(newMesh);
  for (Face f : oldMesh->faces()) {
    outData[fMap[f]] = inData[f];
  }
  for (Face f : oldMesh->boundaryLoops()) {
    outData[fMap[f]] = inData[f];
  }

  return outData;
}

template <class T>
EdgeData<T> HalfedgeMeshDataTransfer::transfer(EdgeData<T>& inData) {

  EdgeData<T> outData(newMesh);
  for (Edge e : oldMesh->edges()) {
    outData[eMap[e]] = inData[e];
  }

  return outData;
}

template <class T>
HalfedgeData<T> HalfedgeMeshDataTransfer::transfer(HalfedgeData<T>& inData) {

  HalfedgeData<T> outData(newMesh);
  for (Halfedge he : oldMesh->allHalfedges()) {
    outData[heMap[he]] = inData[he];
  }

  return outData;
}

template <class T>
CornerData<T> HalfedgeMeshDataTransfer::transfer(CornerData<T>& inData) {

  CornerData<T> outData(newMesh);
  for (Halfedge he : oldMesh->allHalfedges()) {
    outData[heMap[he].corner()] = inData[he.corner()];
  }

  return outData;
}


template <class T>
VertexData<T> HalfedgeMeshDataTransfer::transferBack(VertexData<T>& inData) {

  VertexData<T> outData(oldMesh);
  for (Vertex v : oldMesh->vertices()) {
    outData[v] = inData[vMap[v]];
  }
  return outData;
}


template <class T>
FaceData<T> HalfedgeMeshDataTransfer::transferBack(FaceData<T>& inData) {

  FaceData<T> outData(oldMesh);
  for (Face f : oldMesh->faces()) {
    outData[f] = inData[fMap[f]];
  }
  for (Face f : oldMesh->boundaryLoops()) {
    outData[f] = inData[fMap[f]];
  }
  return outData;
}


template <class T>
EdgeData<T> HalfedgeMeshDataTransfer::transferBack(EdgeData<T>& inData) {

  EdgeData<T> outData(oldMesh);
  for (Edge e : oldMesh->edges()) {
    outData[e] = inData[eMap[e]];
  }
  return outData;
}

template <class T>
HalfedgeData<T> HalfedgeMeshDataTransfer::transferBack(HalfedgeData<T>& inData) {

  HalfedgeData<T> outData(oldMesh);
  for (Halfedge he : oldMesh->allHalfedges()) {
    outData[he] = inData[heMap[he]];
  }
  return outData;
}

template <class T>
CornerData<T> HalfedgeMeshDataTransfer::transferBack(CornerData<T>& inData) {

  CornerData<T> outData(oldMesh);
  for (Halfedge he : oldMesh->allHalfedges()) {
    outData[he.corner()] = inData[heMap[he].corner()];
  }
  return outData;
}


} // namespace geometrycentral
