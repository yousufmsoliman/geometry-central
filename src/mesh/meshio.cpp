#include <iostream>
#include <limits>

#include "geometrycentral/mesh/meshio.h"
#include "geometrycentral/geometry/geometry.h"
#include "geometrycentral/mesh/polygon_soup_mesh.h"

#include "happly.h"

using namespace std;

namespace geometrycentral {

// ======= Input =======

Geometry<Euclidean>* loadMesh_PLY(std::string filename) {

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
  Geometry<Euclidean>* geometry = nullptr;
  HalfedgeMesh* mesh = new HalfedgeMesh(PolygonSoupMesh(faceIndices, vertexPositions), geometry);

  return geometry;
}

Geometry<Euclidean>* loadMesh_OBJ(std::string filename) {

  Geometry<Euclidean>* geometry = nullptr;
  HalfedgeMesh* mesh = new HalfedgeMesh(PolygonSoupMesh(filename), geometry);

  return geometry;
}

Geometry<Euclidean>* loadMesh(std::string filename, std::string type) {

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
    return loadMesh_OBJ(filename);
  } else if (type == "ply") {
    return loadMesh_PLY(filename);
  } else {
    if(typeGiven) {
      throw std::runtime_error("Did not recognize mesh file type " + type);
    } else {
      throw std::runtime_error("Could not detect file type to load mesh from " + filename);
    }
  }

  return nullptr;
}


// ======= Output =======

bool WavefrontOBJ::write(string filename, Geometry<Euclidean>& geometry) {
  ofstream out;
  if (!openStream(out, filename)) return false;

  writeHeader(out, geometry);
  out << "# texture coordinates: NO" << endl;
  cout << endl;

  writeVertices(out, geometry);

  bool useTexCoords = false;
  writeFaces(out, geometry, useTexCoords);

  return true;
}

bool WavefrontOBJ::write(string filename, Geometry<Euclidean>& geometry, CornerData<Vector2>& texcoords) {
  ofstream out;
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

bool WavefrontOBJ::openStream(ofstream& out, string filename) {
  out.open(filename);

  if (!out.is_open()) {
    return false;
  }

  // Use full precision---yes, this makes files bigger, but it also
  // means that saving and then re-loading files doesn't result in
  // unexpected behavior.
  out.precision(numeric_limits<double>::max_digits10);

  return true;
}

void WavefrontOBJ::writeHeader(ofstream& out, Geometry<Euclidean>& geometry) {
  out << "# Mesh exported from GeometryCentral" << endl;
  out << "#  vertices: " << geometry.mesh.nVertices() << endl;
  out << "#     edges: " << geometry.mesh.nEdges() << endl;
  out << "#     faces: " << geometry.mesh.nFaces() << endl;
}

void WavefrontOBJ::writeVertices(ofstream& out, Geometry<Euclidean>& geometry) {
  HalfedgeMesh& mesh(geometry.mesh);

  for (VertexPtr v : mesh.vertices()) {
    Vector3 p = geometry.position(v);
    out << "v " << p.x << " " << p.y << " " << p.z << endl;
  }
}

void WavefrontOBJ::writeTexCoords(ofstream& out, Geometry<Euclidean>& geometry, CornerData<Vector2>& texcoords) {
  HalfedgeMesh& mesh(geometry.mesh);

  for (CornerPtr c : mesh.corners()) {
    if (c.halfedge().isReal()) {
      Vector2 z = texcoords[c];
      out << "vt " << z.x << " " << z.y << endl;
    }
  }
}

void WavefrontOBJ::writeFaces(ofstream& out, Geometry<Euclidean>& geometry, bool useTexCoords) {
  HalfedgeMesh& mesh(geometry.mesh);

  // Get vertex indices
  VertexData<size_t> indices = mesh.getVertexIndices();

  if (useTexCoords) {
    // Get corner indices
    CornerData<size_t> cIndices = mesh.getCornerIndices();

    for (FacePtr f : mesh.faces()) {
      out << "f";
      for (CornerPtr c : f.adjacentCorners()) {
        out << " " << indices[c.vertex()] + 1 << "/" << cIndices[c] + 1; // OBJ uses 1-based indexing
      }
      out << endl;
    }
  } else {
    for (FacePtr f : mesh.faces()) {
      out << "f";
      for (HalfedgePtr h : f.adjacentHalfedges()) {
        out << " " << indices[h.vertex()] + 1; // OBJ uses 1-based indexing
      }
      out << endl;
    }
  }
}

bool PLY::write(std::string filename, Geometry<Euclidean>& geometry, VertexData<Vector3> colors) {
  ofstream out;
  out.open(filename);
  if (!out) {
    cerr << "Could no open file for writing: " << filename << endl;
    return false;
  }
  out.precision(numeric_limits<double>::max_digits10);

  HalfedgeMesh* mesh = geometry.getMesh();
  size_t nVert = mesh->nVertices();
  size_t nFace = mesh->nFaces();

  // Header
  out << "ply" << endl;
  out << "format ascii 1.0" << endl;
  out << "element vertex " << nVert << endl;
  out << "property float x" << endl;
  out << "property float y" << endl;
  out << "property float z" << endl;
  out << "property uchar red" << endl;
  out << "property uchar green" << endl;
  out << "property uchar blue" << endl;
  out << "element face " << nFace << endl;
  out << "property list uchar int vertex_index" << endl;
  out << "end_header" << endl;

  auto round255 = [&](double v) { return static_cast<int>(clamp(std::round(v * 255), 0.0, 255.)); };

  // Vertices
  for (VertexPtr v : mesh->vertices()) {
    out << geometry.position(v).x << " ";
    out << geometry.position(v).y << " ";
    out << geometry.position(v).z << " ";
    out << round255(colors[v].x) << " ";
    out << round255(colors[v].y) << " ";
    out << round255(colors[v].z) << endl;
  }

  // Faces
  VertexData<size_t> vInd = mesh->getVertexIndices();
  for (FacePtr f : mesh->faces()) {
    out << f.degree();
    for (VertexPtr v : f.adjacentVertices()) {
      out << " " << vInd[v];
    }
    out << endl;
  }

  out.close();

  return true;
}

} // namespace geometrycentral
