#pragma once

// Transfer data between halfedge meshes after a mesh.copy()
// (so we assume exactly discrete correspondence)
// Don't make this yourself; the maps get populated when you call
// HalfedgeMesh::copy();

#include <unordered_map>

namespace geometrycentral {

class HalfedgeMeshDataTransfer {

public:
  HalfedgeMeshDataTransfer();
  HalfedgeMeshDataTransfer(HalfedgeMesh* oldMesh, HalfedgeMesh* newMesh);

  // Transfer data from old to new
  template <class T>
  VertexData<T> transfer(VertexData<T>& inData);
  template <class T>
  FaceData<T> transfer(FaceData<T>& inData);
  template <class T>
  EdgeData<T> transfer(EdgeData<T>& inData);
  template <class T>
  HalfedgeData<T> transfer(HalfedgeData<T>& inData);
  template <class T>
  CornerData<T> transfer(CornerData<T>& inData);

  // Transfer data from new to old
  template <class T>
  VertexData<T> transferBack(VertexData<T>& inData);
  template <class T>
  FaceData<T> transferBack(FaceData<T>& inData);
  template <class T>
  EdgeData<T> transferBack(EdgeData<T>& inData);
  template <class T>
  HalfedgeData<T> transferBack(HalfedgeData<T>& inData);
  template <class T>
  CornerData<T> transferBack(CornerData<T>& inData);

  // Maps to transfer pointers map[old] --> new
  HalfedgeData<Halfedge> heMap;
  VertexData<Vertex> vMap;
  EdgeData<Edge> eMap;
  FaceData<Face> fMap;


  // Maps to transfer pointers map[new] --> old
  HalfedgeData<Halfedge> heMapBack;
  VertexData<Vertex> vMapBack;
  EdgeData<Edge> eMapBack;
  FaceData<Face> fMapBack;
  void generateReverseMaps();

  // Relevant meshes
  HalfedgeMesh* oldMesh = nullptr;
  HalfedgeMesh* newMesh = nullptr;
};

} // namespace geometrycentral
