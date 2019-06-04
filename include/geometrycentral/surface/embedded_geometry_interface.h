#pragma once

#include "geometrycentral/surface/extrinsic_geometry_interface.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/utilities/vector3.h"

#include <Eigen/SparseCore>

namespace geometrycentral {
namespace surface {


class EmbeddedGeometryInterface : public ExtrinsicGeometryInterface {

protected:
  // Constructor is protected, because this class is an interface which is not meant to be instantiated directly.
  // Instantiate it via some realization which encapsulates input data, like VertexPositionGeometry.
  EmbeddedGeometryInterface(HalfedgeMesh& mesh_);
  virtual ~EmbeddedGeometryInterface() {}

public:
  // == Quantities

  // Vertex positions
  inline void requireVertexPositions();
  inline void unrequireVertexPositions();
  VertexData<Vector3> vertexPositions;

protected:
  // Implmentations of quantities from base classes
  virtual void computeEdgeLengths() override;

  DependentQuantityD<VertexData<Vector3>> vertexPositionsQ;
  virtual void computeVertexPositions() = 0;
};


} // namespace surface
} // namespace geometrycentral
