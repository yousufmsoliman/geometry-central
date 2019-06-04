#include "geometrycentral/surface/embedded_geometry_interface.h"

#include <limits>

namespace geometrycentral {
namespace surface {

// clang-format off
EmbeddedGeometryInterface::EmbeddedGeometryInterface(HalfedgeMesh& mesh_) : 
  ExtrinsicGeometryInterface(mesh_),

  vertexPositionsQ      (&vertexPositions,      std::bind(&EmbeddedGeometryInterface::computeVertexPositions, this),       quantities)
  //faceAreasQ                (&faceAreas,                std::bind(&EmbeddedGeometryInterface::computeFaceAreas, this),          {})
  
  {}
// clang-format on

void EmbeddedGeometryInterface::computeEdgeLengths() {
  vertexPositionsQ.ensureHave();

  edgeLengths = EdgeData<double>(mesh);
  for (Edge e : mesh.edges()) {
    edgeLengths[e] = norm(vertexPositions[e.halfedge().vertex()] - vertexPositions[e.halfedge().twin().vertex()]);
  }
}


void EmbeddedGeometryInterface::requireVertexPositions() { vertexPositionsQ.require(); }
void EmbeddedGeometryInterface::unrequireVertexPositions() { vertexPositionsQ.unrequire(); }

} // namespace surface
} // namespace geometrycentral
