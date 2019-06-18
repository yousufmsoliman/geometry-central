#include "geometrycentral/surface/embedded_geometry_interface.h"

#include <limits>

using std::cout;
using std::endl;

namespace geometrycentral {
namespace surface {

// clang-format off
EmbeddedGeometryInterface::EmbeddedGeometryInterface(HalfedgeMesh& mesh_) : 
  ExtrinsicGeometryInterface(mesh_),

  vertexPositionsQ      (&vertexPositions,      std::bind(&EmbeddedGeometryInterface::computeVertexPositions, this),        quantities),
  faceNormalsQ          (&faceNormals,          std::bind(&EmbeddedGeometryInterface::computeFaceNormals, this),            quantities),
  vertexNormalsQ        (&vertexNormals,        std::bind(&EmbeddedGeometryInterface::computeVertexNormals, this),          quantities),
  faceTangentBasisQ     (&faceTangentBasis,     std::bind(&EmbeddedGeometryInterface::computeFaceTangentBasis, this),       quantities),
  vertexTangentBasisQ   (&vertexTangentBasis,   std::bind(&EmbeddedGeometryInterface::computeVertexTangentBasis, this),     quantities)
  
  {}
// clang-format on


// Edge lengths
double EmbeddedGeometryInterface::edgeLength(Edge e) const {
  return norm(vertexPosition(e.halfedge().vertex()) - vertexPosition(e.halfedge().twin().vertex()));
}
void EmbeddedGeometryInterface::computeEdgeLengths() {
  vertexPositionsQ.ensureHave();

  edgeLengths = EdgeData<double>(mesh);
  for (Edge e : mesh.edges()) {
    edgeLengths[e] = norm(vertexPositions[e.halfedge().vertex()] - vertexPositions[e.halfedge().twin().vertex()]);
  }
}

// Edge dihedral angles
void EmbeddedGeometryInterface::computeEdgeDihedralAngles() {
  vertexPositionsQ.ensureHave();
  faceNormalsQ.ensureHave();

  edgeDihedralAngles = EdgeData<double>(mesh, 0.);
  for (Edge e : mesh.edges()) {
    if (e.isBoundary()) continue;

    Vector3 N1 = faceNormals[e.halfedge().face()];
    Vector3 N2 = faceNormals[e.halfedge().twin().face()];
    Vector3 pTail = vertexPositions[e.halfedge().vertex()];
    Vector3 pTip = vertexPositions[e.halfedge().twin().vertex()];
    Vector3 edgeDir = unit(-pTip + pTail);

    edgeDihedralAngles[e] = atan2(dot(edgeDir, cross(N1, N2)), dot(N1, N2));
  }
}
void EmbeddedGeometryInterface::requireVertexPositions() { vertexPositionsQ.require(); }
void EmbeddedGeometryInterface::unrequireVertexPositions() { vertexPositionsQ.unrequire(); }


// Face normal
Vector3 EmbeddedGeometryInterface::faceNormal(Face f) const {
  // For general polygons, take the sum of the cross products at each corner
  Vector3 normalSum = Vector3::zero();
  for (Halfedge heF : f.adjacentHalfedges()) {

    // Gather vertex positions for next three vertices
    Halfedge he = heF;
    Vector3 pA = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pB = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = vertexPositions[he.vertex()];

    normalSum += cross(pB - pA, pC - pA);

    // In the special case of a triangle, there is no need to to repeat at all three corners; the result will be the
    // same
    if (he.next() == heF) break;
  }

  Vector3 normal = unit(normalSum);
  return normal;
}
void EmbeddedGeometryInterface::computeFaceNormals() {
  vertexPositionsQ.ensureHave();

  faceNormals = FaceData<Vector3>(mesh);

  for (Face f : mesh.faces()) {

    // For general polygons, take the sum of the cross products at each corner
    Vector3 normalSum = Vector3::zero();
    for (Halfedge heF : f.adjacentHalfedges()) {

      // Gather vertex positions for next three vertices
      Halfedge he = heF;
      Vector3 pA = vertexPositions[he.vertex()];
      he = he.next();
      Vector3 pB = vertexPositions[he.vertex()];
      he = he.next();
      Vector3 pC = vertexPositions[he.vertex()];

      normalSum += cross(pB - pA, pC - pA);

      // In the special case of a triangle, there is no need to to repeat at all three corners; the result will be the
      // same
      if (he.next() == heF) break;
    }

    Vector3 normal = unit(normalSum);
    faceNormals[f] = normal;
  }
}
void EmbeddedGeometryInterface::requireFaceNormals() { faceNormalsQ.require(); }
void EmbeddedGeometryInterface::unrequireFaceNormals() { faceNormalsQ.unrequire(); }

// Vertex normal
void EmbeddedGeometryInterface::computeVertexNormals() {
  faceNormalsQ.ensureHave();
  cornerAnglesQ.ensureHave();

  vertexNormals = VertexData<Vector3>(mesh);

  for (Vertex v : mesh.vertices()) {
    Vector3 normalSum = Vector3::zero();

    for (Corner c : v.adjacentCorners()) {
      Vector3 normal = faceNormals[c.face()];
      double weight = cornerAngles[c];

      normalSum += weight * normal;
    }

    vertexNormals[v] = unit(normalSum);
  }
}
void EmbeddedGeometryInterface::requireVertexNormals() { vertexNormalsQ.require(); }
void EmbeddedGeometryInterface::unrequireVertexNormals() { vertexNormalsQ.unrequire(); }

// Face tangent basis
void EmbeddedGeometryInterface::computeFaceTangentBasis() {
  vertexPositionsQ.ensureHave();
  faceNormalsQ.ensureHave();
  halfedgeVectorsInFaceQ.ensureHave();

  faceTangentBasis = FaceData<std::array<Vector3, 2>>(mesh);

  for (Face f : mesh.faces()) {

    // For general polygons, take the average of each edge vector projected to tangent plane
    Vector3 basisXSum = Vector3::zero();
    Vector3 N = faceNormals[f];
    bool isTriangular = f.isTriangle();
    for (Halfedge heF : f.adjacentHalfedges()) {

      Vector3 eVec = vertexPositions[heF.twin().vertex()] - vertexPositions[heF.vertex()];
      eVec.removeComponent(N);

      // TODO can surely do this with less trig
      double angle = halfedgeVectorsInFace[heF].arg();
      Vector3 eVecX = eVec.rotate_around(N, -angle);

      basisXSum += eVecX;

      // In the special case of a triangle, there is no need to to repeat at all three corners; the result will be the
      // same
      if (isTriangular) break;
    }

    Vector3 basisX = unit(basisXSum);
    Vector3 basisY = cross(N, basisX);
    faceTangentBasis[f] = {{basisX, basisY}};
  }
}
void EmbeddedGeometryInterface::requireFaceTangentBasis() { faceTangentBasisQ.require(); }
void EmbeddedGeometryInterface::unrequireFaceTangentBasis() { faceTangentBasisQ.unrequire(); }

// Vertex tangent basis
void EmbeddedGeometryInterface::computeVertexTangentBasis() {
  vertexPositionsQ.ensureHave();
  vertexNormalsQ.ensureHave();
  halfedgeVectorsInVertexQ.ensureHave();

  vertexTangentBasis = VertexData<std::array<Vector3, 2>>(mesh);

  for (Vertex v : mesh.vertices()) {

    // For general polygons, take the average of each edge vector projected to tangent plane
    Vector3 basisXSum = Vector3::zero();
    Vector3 N = vertexNormals[v];
    for (Halfedge he : v.outgoingHalfedges()) {

      Vector3 eVec = vertexPositions[he.twin().vertex()] - vertexPositions[he.vertex()];
      eVec.removeComponent(N);

      // TODO can surely do this with less trig
      double angle = halfedgeVectorsInVertex[he].arg();
      Vector3 eVecX = eVec.rotate_around(N, -angle);

      basisXSum += eVecX;
    }

    Vector3 basisX = unit(basisXSum);
    Vector3 basisY = cross(N, basisX);
    vertexTangentBasis[v] = {{basisX, basisY}};
  }
}
void EmbeddedGeometryInterface::requireVertexTangentBasis() { vertexTangentBasisQ.require(); }
void EmbeddedGeometryInterface::unrequireVertexTangentBasis() { vertexTangentBasisQ.unrequire(); }

} // namespace surface
} // namespace geometrycentral
