#include "geometrycentral/surface/geometry_base.h"

namespace geometrycentral {
namespace surface {


// clang-format off
GeometryBase::GeometryBase(HalfedgeMesh& mesh_)
    : mesh(mesh_),
      
  // Construct the dependency graph of managed quantities and their callbacks

  vertexIndicesQ           (&vertexIndices,         std::bind(&GeometryBase::computeVertexIndices, this),           {}),
  interiorVertexIndicesQ   (&interiorVertexIndices, std::bind(&GeometryBase::computeInteriorVertexIndices, this),   {}),
  edgeIndicesQ             (&edgeIndices,           std::bind(&GeometryBase::computeEdgeIndices, this),             {}),
  halfedgeIndicesQ         (&halfedgeIndices,       std::bind(&GeometryBase::computeHalfedgeIndices, this),         {}),
  cornerIndicesQ           (&cornerIndices,         std::bind(&GeometryBase::computeCornerIndices, this),           {}),
  faceIndicesQ             (&faceIndices,           std::bind(&GeometryBase::computeFaceIndices, this),             {})
  //boundaryLoopIndicesQ     (&boundaryLoopIndices,   std::bind(&GeometryBase::computeBoundaryLoopIndices, this),     {})

  {

    quantities.push_back(&vertexIndicesQ);

    std::cout << "constructed geometry_base!" << std::endl;
  }
// clang-format on

GeometryBase::~GeometryBase() {}

void GeometryBase::refreshQuantities() {
  for (DependentQuantity* q : quantities) {
    q->computed = false;
  }
  for (DependentQuantity* q : quantities) {
    q->ensureHaveIfRequired();
  }
}

void GeometryBase::purgeQuantities() {
  for (DependentQuantity* q : quantities) {
    q->clearIfNotRequired();
  }
}


// == Indices

void GeometryBase::computeVertexIndices() { vertexIndices = mesh.getVertexIndices(); }
void GeometryBase::computeInteriorVertexIndices() { interiorVertexIndices = mesh.getInteriorVertexIndices(); }
void GeometryBase::computeEdgeIndices() { edgeIndices = mesh.getEdgeIndices(); }
void GeometryBase::computeHalfedgeIndices() { halfedgeIndices = mesh.getHalfedgeIndices(); }
void GeometryBase::computeCornerIndices() { cornerIndices = mesh.getCornerIndices(); }
void GeometryBase::computeFaceIndices() { faceIndices = mesh.getFaceIndices(); }


} // namespace surface
} // namespace geometrycentral
