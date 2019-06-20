#include "geometrycentral/surface/meshio.h"

#include "geometrycentral/surface/halfedge_containers.h"
#include "geometrycentral/surface/halfedge_factories.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/polygon_soup_mesh.h"

#include "happly.h"

#include <iostream>
#include <limits>

using std::cout;
using std::endl;

namespace geometrycentral {
namespace surface {

// ======= Input =======

// Mesh loader helpers
namespace {
std::tuple<std::unique_ptr<HalfedgeMesh>, std::unique_ptr<VertexPositionGeometry>> loadMesh_PLY(std::string filename,
                                                                                                bool verbose) {

  happly::PLYData plyData(filename);

  // === Get vertex positions
  std::vector<std::array<double, 3>> rawPos = plyData.getVertexPositions();
  std::vector<Vector3> vertexPositions(rawPos.size());
  for (size_t i = 0; i < rawPos.size(); i++) {
    vertexPositions[i][0] = rawPos[i][0];
    vertexPositions[i][1] = rawPos[i][1];
    vertexPositions[i][2] = rawPos[i][2];
  }

  // Get face list
  std::vector<std::vector<size_t>> faceIndices = plyData.getFaceIndices();

  // === Build the mesh objects
  return makeHalfedgeAndGeometry(faceIndices, vertexPositions, verbose);
}

std::tuple<std::unique_ptr<HalfedgeMesh>, std::unique_ptr<VertexPositionGeometry>> loadMesh_OBJ(std::string filename,
                                                                                                bool verbose) {
  PolygonSoupMesh soup(filename);
  return makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates, verbose);
}
} // namespace

std::tuple<std::unique_ptr<HalfedgeMesh>, std::unique_ptr<VertexPositionGeometry>>
loadMesh(std::string filename, bool verbose, std::string type) {

  // Check if file exists
  std::ifstream testStream(filename);
  if (!testStream) {
    throw std::runtime_error("Could not load mesh; file does not exist: " + filename);
  }
  testStream.close();

  // Attempt to detect filename
  bool typeGiven = type != "";
  std::string::size_type sepInd = filename.rfind('.');
  if (!typeGiven) {
    if (sepInd != std::string::npos) {
      std::string extension;
      extension = filename.substr(sepInd + 1);

      // Convert to all lowercase
      std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
      type = extension;
    }
  }


  if (type == "obj") {
    return loadMesh_OBJ(filename, verbose);
  } else if (type == "ply") {
    return loadMesh_PLY(filename, verbose);
  } else {
    if (typeGiven) {
      throw std::runtime_error("Did not recognize mesh file type " + type);
    } else {
      throw std::runtime_error("Could not detect file type to load mesh from " + filename);
    }
  }
}


// connectivity loader helpers
namespace {
std::unique_ptr<HalfedgeMesh> loadConnectivity_PLY(std::string filename, bool verbose) {

  happly::PLYData plyData(filename);

  // Get face list
  std::vector<std::vector<size_t>> faceIndices = plyData.getFaceIndices();

  // === Build the mesh objects
  return std::unique_ptr<HalfedgeMesh>(new HalfedgeMesh(faceIndices, verbose));
}

std::unique_ptr<HalfedgeMesh> loadConnectivity_OBJ(std::string filename, bool verbose) {
  // TODO this will fail unless the obj file has vertex listings, which is not strictly needed to load connectivity. I'm
  // not sure if that's really a valid .obj file, but nonetheless this function could certainly process such .obj files.
  PolygonSoupMesh soup(filename);
  return std::unique_ptr<HalfedgeMesh>(new HalfedgeMesh(soup.polygons, verbose));
}
} // namespace


std::unique_ptr<HalfedgeMesh> loadConnectivity(std::string filename, bool verbose, std::string type) {

  // Check if file exists
  std::ifstream testStream(filename);
  if (!testStream) {
    throw std::runtime_error("Could not load mesh; file does not exist: " + filename);
  }
  testStream.close();

  // Attempt to detect filename
  bool typeGiven = type != "";
  std::string::size_type sepInd = filename.rfind('.');
  if (!typeGiven) {
    if (sepInd != std::string::npos) {
      std::string extension;
      extension = filename.substr(sepInd + 1);

      // Convert to all lowercase
      std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
      type = extension;
    }
  }


  if (type == "obj") {
    return loadConnectivity_OBJ(filename, verbose);
  } else if (type == "ply") {
    return loadConnectivity_PLY(filename, verbose);
  } else {
    if (typeGiven) {
      throw std::runtime_error("Did not recognize mesh file type " + type);
    } else {
      throw std::runtime_error("Could not detect file type to load mesh from " + filename);
    }
  }
}


// ======= Output =======

/*
bool WavefrontOBJ::write(std::string filename, GeometryEuclidean>& geometry) {
  std::ofstream out;
  if (!openStream(out, filename)) return false;

  writeHeader(out, geometry);
  out << "# texture coordinates: NO" << endl;
  cout << endl;

  writeVertices(out, geometry);

  bool useTexCoords = false;
  writeFaces(out, geometry, useTexCoords);

  return true;
}

bool WavefrontOBJ::write(std::string filename, Geometry<Euclidean>& geometry, CornerData<Vector2>& texcoords) {
  std::ofstream out;
  if (!openStream(out, filename)) return false;

  writeHeader(out, geometry);
  out << "# texture coordinates: YES" << endl;
  cout << endl;

  writeVertices(out, geometry);
  writeTexCoords(out, geometry, texcoords);

  bool useTexCoords = true;
  writeFaces(out, geometry, useTexCoords);

  return true;
}

bool WavefrontOBJ::openStream(std::ofstream& out, std::string filename) {
  out.open(filename);

  if (!out.is_open()) {
    return false;
  }

  // Use full precision---yes, this makes files bigger, but it also
  // means that saving and then re-loading files doesn't result in
  // unexpected behavior.
  out.precision(std::numeric_limits<double>::max_digits10);

  return true;
}

void WavefrontOBJ::writeHeader(std::ofstream& out, Geometry<Euclidean>& geometry) {
  out << "# Mesh exported from GeometryCentral" << endl;
  out << "#  vertices: " << geometry.mesh.nVertices() << endl;
  out << "#     edges: " << geometry.mesh.nEdges() << endl;
  out << "#     faces: " << geometry.mesh.nFaces() << endl;
}

void WavefrontOBJ::writeVertices(std::ofstream& out, Geometry<Euclidean>& geometry) {
  HalfedgeMesh& mesh(geometry.mesh);

  for (Vertex v : mesh.vertices()) {
    Vector3 p = geometry.position(v);
    out << "v " << p.x << " " << p.y << " " << p.z << endl;
  }
}

void WavefrontOBJ::writeTexCoords(std::ofstream& out, Geometry<Euclidean>& geometry, CornerData<Vector2>& texcoords) {
  HalfedgeMesh& mesh(geometry.mesh);

  for (Corner c : mesh.corners()) {
    Vector2 z = texcoords[c];
    out << "vt " << z.x << " " << z.y << endl;
  }
}

void WavefrontOBJ::writeFaces(std::ofstream& out, Geometry<Euclidean>& geometry, bool useTexCoords) {
  HalfedgeMesh& mesh(geometry.mesh);

  // Get vertex indices
  VertexData<size_t> indices = mesh.getVertexIndices();

  if (useTexCoords) {
    // Get corner indices
    CornerData<size_t> cIndices = mesh.getCornerIndices();

    for (Face f : mesh.faces()) {
      out << "f";
      for (Corner c : f.adjacentCorners()) {
        out << " " << indices[c.vertex()] + 1 << "/" << cIndices[c] + 1; // OBJ uses 1-based indexing
      }
      out << endl;
    }
  } else {
    for (Face f : mesh.faces()) {
      out << "f";
      for (Halfedge h : f.adjacentHalfedges()) {
        out << " " << indices[h.vertex()] + 1; // OBJ uses 1-based indexing
      }
      out << endl;
    }
  }
}
*/


} // namespace surface
} // namespace geometrycentral
