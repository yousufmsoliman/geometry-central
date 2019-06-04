#include "geometrycentral/surface/geometry_base.h"

namespace geometrycentral {
namespace surface {


// clang-format off
GeometryBase::GeometryBase(HalfedgeMesh& mesh_, std::vector<DependentQuantity*> childQuantities)
    : mesh(mesh_),
      
  // Construct the dependency graph of managed quantities and their callbacks

  vertexIndicesQ           (&vertexIndices,         std::bind(&GeometryBase::computeVertexIndices, this),           {}),
  interiorVertexIndicesQ   (&interiorVertexIndices, std::bind(&GeometryBase::computeInteriorVertexIndices, this),   {}),
  edgeIndicesQ             (&edgeIndices,           std::bind(&GeometryBase::computeEdgeIndices, this),             {}),
  halfedgeIndicesQ         (&halfedgeIndices,       std::bind(&GeometryBase::computeHalfedgeIndices, this),         {}),
  cornerIndicesQ           (&cornerIndices,         std::bind(&GeometryBase::computeCornerIndices, this),           {}),
  faceIndicesQ             (&faceIndices,           std::bind(&GeometryBase::computeFaceIndices, this),             {}),
  boundaryLoopIndicesQ     (&boundaryLoopIndices,   std::bind(&GeometryBase::computeBoundaryLoopIndices, this),     {})

  {
    quantities.push_back(&vertexIndicesQ);
    quantities.push_back(&interiorVertexIndicesQ);
    quantities.push_back(&edgeIndicesQ);
    quantities.push_back(&halfedgeIndicesQ);
    quantities.push_back(&cornerIndicesQ);
    quantities.push_back(&faceIndicesQ);
    quantities.push_back(&boundaryLoopIndicesQ);
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
  
std::unique_ptr<GeometryBase> GeometryBase::reinterpretTo(HalfedgeMesh& targetMesh) {
  std::unique_ptr<GeometryBase> newGeom(new GeometryBase(targetMesh));
  return newGeom; 
}


// == Indices

void GeometryBase::computeVertexIndices() { vertexIndices = mesh.getVertexIndices(); }
void GeometryBase::computeInteriorVertexIndices() { interiorVertexIndices = mesh.getInteriorVertexIndices(); }
void GeometryBase::computeEdgeIndices() { edgeIndices = mesh.getEdgeIndices(); }
void GeometryBase::computeHalfedgeIndices() { halfedgeIndices = mesh.getHalfedgeIndices(); }
void GeometryBase::computeCornerIndices() { cornerIndices = mesh.getCornerIndices(); }
void GeometryBase::computeFaceIndices() { faceIndices = mesh.getFaceIndices(); }
void GeometryBase::computeBoundaryLoopIndices() { boundaryLoopIndices = mesh.getBoundaryLoopIndices(); }


} // namespace surface
} // namespace geometrycentral
