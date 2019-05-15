#include "geometrycentral/mesh/ply_wrapper.h"

#include "geometrycentral/mesh/polygon_soup_mesh.h"

#include <cstring>
#include <fstream>
#include <iostream>


using std::cerr;
using std::cout;
using std::endl;

namespace geometrycentral {


PlyHalfedgeMeshData::PlyHalfedgeMeshData(std::string filename_, bool verbose_)
    : filename(filename_), verbose(verbose_) {

  plyData = new happly::PLYData(filename, verbose);
}

PlyHalfedgeMeshData::PlyHalfedgeMeshData(Geometry<Euclidean>* geometry_)
    : mesh(geometry_->getMesh()), geometry(geometry_) {

  plyData = new happly::PLYData();

  // Vertex positions. Need to convert type for happly.
  std::vector<Vector3> vertexPositions = geometry->getVertexPositionList();
  std::vector<std::array<double, 3>> vertexPositionsA(vertexPositions.size());
  for (size_t v = 0; v < vertexPositions.size(); v++) {
    for (size_t i = 0; i < 3; i++) {
      vertexPositionsA[v][i] = vertexPositions[v][i];
    }
  }
  plyData->addVertexPositions(vertexPositionsA);

  // Face indices
  std::vector<std::vector<size_t>> faceIndices = mesh->getPolygonSoupFaces();
  plyData->addFaceIndices(faceIndices);
}

PlyHalfedgeMeshData::~PlyHalfedgeMeshData() { safeDelete(plyData); }


Geometry<Euclidean>* PlyHalfedgeMeshData::getMesh() {

  // Return the cached mesh if this has already been called
  if (mesh != nullptr) {
    return geometry;
  }

  // === Get vertex positions
  std::vector<std::array<double, 3>> rawPos = plyData->getVertexPositions(vertexName);
  std::vector<Vector3> vertexPositions(rawPos.size());
  for (size_t i = 0; i < rawPos.size(); i++) {
    vertexPositions[i][0] = rawPos[i][0];
    vertexPositions[i][1] = rawPos[i][1];
    vertexPositions[i][2] = rawPos[i][2];
  }


  // Get face list
  std::vector<std::vector<size_t>> faceIndices = plyData->getFaceIndices();

  // === Build the mesh objects
  mesh = new HalfedgeMesh(PolygonSoupMesh(faceIndices, vertexPositions), geometry);

  return geometry;
}

VertexData<Vector3> PlyHalfedgeMeshData::getVertexColors() {

  VertexData<Vector3> color(mesh);

  try {
    // Try uchar first
    VertexData<unsigned char> r = getVertexProperty<unsigned char>("red");
    VertexData<unsigned char> g = getVertexProperty<unsigned char>("green");
    VertexData<unsigned char> b = getVertexProperty<unsigned char>("blue");
    for (Vertex v : mesh->vertices()) {
      color[v][0] = r[v] / 255.0;
      color[v][1] = g[v] / 255.0;
      color[v][2] = b[v] / 255.0;
    }
    return color;


  } catch (std::runtime_error orig_e) {

    // If that doesn't work, try float
    try {
      VertexData<double> r = getVertexProperty<double>("red");
      VertexData<double> g = getVertexProperty<double>("green");
      VertexData<double> b = getVertexProperty<double>("blue");
      for (Vertex v : mesh->vertices()) {
        color[v][0] = r[v];
        color[v][1] = g[v];
        color[v][2] = b[v];
      }
      return color;
    } catch (std::runtime_error second_e) {

      // Initial error will be less confusing, since most expect
      throw std::runtime_error("Could not find vertex colors in PLY file, as uchar or float");
    }
  }
}

void PlyHalfedgeMeshData::write(std::string filename) { plyData->write(filename, happly::DataFormat::Binary); }
} // namespace geometrycentral
