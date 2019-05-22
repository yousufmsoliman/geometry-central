#pragma once

#include "geometrycentral/surface/geometry.h"

#include "happly.h"

#include <fstream>
#include <iostream>
#include <string>

// === Reader and writing classes supporting direct interop between the GeometryCentral mesh types and the .ply file
// === format.

namespace geometrycentral {
namespace surface {


class PlyHalfedgeMeshData {

public:

  // Construct by reading from file, mapping the elements on to an existing mesh.
  // To simultaneously read the mesh encoded by the file, see the static method loadMeshAndData below.
  PlyHalfedgeMeshData(HalfedgeMesh& mesh_, std::string filename, bool verbose = false);
 
  // Construct a data object. Connectivity will be written automatically, and other data fields can be added.
  PlyHalfedgeMeshData(HalfedgeMesh& mesh_, bool verbose = false);
  
  // Convenience factory method to simultaneously read the mesh from a file and
  static std::tuple<std::unique_ptr<HalfedgeMesh>, std::unique_ptr<PlyHalfedgeMeshData>>
  loadMeshAndData(std::string filename, bool verbose = false);

  // The mesh on which the properties in this file are presumed to exist
  HalfedgeMesh& mesh; 



  // === Get properties as geometrycentral types.
  // Will fail with an error if not possible.
  
  // Registers vertex posititions from a geometry object with the file
  Geometry<Euclidean>& geometry getGeometry();
  
  template <class T>
  VertexData<T> getVertexProperty(std::string propertyName);

  template <class T>
  FaceData<T> getFaceProperty(std::string propertyName);

  // === Convenience getters

  // Looks for vertex colors, either as a uchar or a float
  VertexData<Vector3> getVertexColors();

  // === Set properties as geometrycentral types.
 
  // Registers vertex posititions from a geometry object with the file
  void addGeometry(Geometry<Euclidean>& geometry);

  template <class T>
  void addVertexProperty(std::string propertyName, VertexData<T>& vData);

  template <class T>
  void addFaceProperty(std::string propertyName, FaceData<T>& fData);

  // Write this object out to file
  void write(std::string filename);

  // The underlying reader/writer object
  happly::PLYData& plyData;

private:

  // File data

  bool isBinary;
  float version;

  // Options
  bool verbose;
  std::string vertexName = "vertex";
  std::string faceName = "face";
  std::string edgeName = "edge";
  std::string halfedgeName = "halfedge";
  std::string boundaryLoopName = "boundaryloop";
};

// === Implementations
template <class T>
VertexData<T> PlyHalfedgeMeshData::getVertexProperty(std::string propertyName) {
  if (mesh == nullptr) {
    getMesh();
  }

  std::vector<T> rawData = plyData->getElement(vertexName).getProperty<T>(propertyName);

  if (rawData.size() != mesh->nVertices()) {
    throw std::runtime_error("Property " + propertyName + " does not have size equal to number of vertices");
  }

  VertexData<T> result(mesh);
  for (size_t i = 0; i < mesh->nVertices(); i++) {
    result[mesh->vertex(i)] = rawData[i];
  }

  return result;
}


template <class T>
FaceData<T> PlyHalfedgeMeshData::getFaceProperty(std::string propertyName) {
  if (mesh == nullptr) {
    getMesh();
  }

  std::vector<T> rawData = plyData->getElement(faceName).getProperty<T>(propertyName);

  if (rawData.size() != mesh->nVertices()) {
    throw std::runtime_error("Property " + propertyName + " does not have size equal to number of vertices");
  }

  VertexData<T> result(mesh);
  for (size_t i = 0; i < mesh->nVertices(); i++) {
    result[mesh->vertex(i)] = rawData[i];
  }

  return result;
}

template <class T>
void PlyHalfedgeMeshData::addVertexProperty(std::string propertyName, VertexData<T>& vData) {

  std::vector<T> vec;
  for (Vertex v : mesh->vertices()) {
    vec.push_back(vData[v]);
  }

  plyData->getElement(vertexName).addProperty<T>(propertyName, vec);
}


template <class T>
void PlyHalfedgeMeshData::addFaceProperty(std::string propertyName, FaceData<T>& fData) {

  std::vector<T> vec;
  for (Face f : mesh->faces()) {
    vec.push_back(fData[f]);
  }

  plyData->getElement(faceName).addProperty<T>(propertyName, vec);
}

} // namespace surface
} // namespace geometrycentral
