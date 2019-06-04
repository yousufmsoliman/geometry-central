#pragma once

#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/utilities/dependent_quantity.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

namespace geometrycentral {
namespace surface {

class GeometryBase {

public:
  GeometryBase(HalfedgeMesh& mesh, std::vector<DependentQuantity*> childQuantities = {});
  virtual ~GeometryBase();

  // == Members
  HalfedgeMesh& mesh;


  // == Utility methods

  // Recompute all require'd quantities from input data. Call this after e.g. repositioning a vertex or mutating the
  // mesh
  void refreshQuantities();

  // Clear out any cached quantities which were previously computed but are not currently required.
  void purgeQuantities();

  // Construct a geometry object on another mesh identical to this one
  // TODO move this to exist in realizations only
  std::unique_ptr<GeometryBase> reinterpretTo(HalfedgeMesh& targetMesh);


  // === Quantities


  // == Indices

  VertexData<size_t> vertexIndices;
  inline void requireVertexIndices() { vertexIndicesQ.require(); }
  inline void unrequireVertexIndices() { vertexIndicesQ.unrequire(); }

  VertexData<size_t> interiorVertexIndices;
  inline void requireInteriorVertexIndices() { interiorVertexIndicesQ.require(); }
  inline void unrequireInteriorVertexIndices() { interiorVertexIndicesQ.unrequire(); }

  EdgeData<size_t> edgeIndices;
  inline void requireEdgeIndices() { edgeIndicesQ.require(); }
  inline void unrequireEdgeIndices() { edgeIndicesQ.unrequire(); }

  HalfedgeData<size_t> halfedgeIndices;
  inline void requireHalfedgeIndices() { halfedgeIndicesQ.require(); }
  inline void unrequireHalfedgeIndices() { halfedgeIndicesQ.unrequire(); }

  CornerData<size_t> cornerIndices;
  inline void requireCornerIndices() { cornerIndicesQ.require(); }
  inline void unrequireCornerIndices() { cornerIndicesQ.unrequire(); }

  FaceData<size_t> faceIndices;
  inline void requireFaceIndices() { faceIndicesQ.require(); }
  inline void unrequireFaceIndices() { faceIndicesQ.unrequire(); }

  BoundaryLoopData<size_t> boundaryLoopIndices;
  inline void requireBoundaryLoopIndices() { boundaryLoopIndicesQ.require(); }
  inline void unrequireBoundaryLoopIndices() { boundaryLoopIndicesQ.unrequire(); }

protected:

  // All of the quantities available (subclasses will also add quantities to this list)
  std::vector<DependentQuantity*> quantities;

  // === Implementation details for quantities

  // == Indices

  DependentQuantityD<VertexData<size_t>> vertexIndicesQ;
  virtual void computeVertexIndices();

  DependentQuantityD<VertexData<size_t>> interiorVertexIndicesQ;
  virtual void computeInteriorVertexIndices();

  DependentQuantityD<EdgeData<size_t>> edgeIndicesQ;
  virtual void computeEdgeIndices();

  DependentQuantityD<HalfedgeData<size_t>> halfedgeIndicesQ;
  virtual void computeHalfedgeIndices();

  DependentQuantityD<CornerData<size_t>> cornerIndicesQ;
  virtual void computeCornerIndices();

  DependentQuantityD<FaceData<size_t>> faceIndicesQ;
  virtual void computeFaceIndices();

  DependentQuantityD<BoundaryLoopData<size_t>> boundaryLoopIndicesQ;
  virtual void computeBoundaryLoopIndices();
};

} // namespace surface
} // namespace geometrycentral
