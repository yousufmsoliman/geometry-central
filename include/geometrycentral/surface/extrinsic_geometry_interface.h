#pragma once

#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/utilities/vector2.h"

#include <Eigen/SparseCore>

namespace geometrycentral {
namespace surface {


class ExtrinsicGeometryInterface : public IntrinsicGeometryInterface {

protected:
  // Constructor is protected, because this class is an interface which is not meant to be instantiated directly.
  // Instantiate it via some realization which encapsulates input data, like VertexPositionGeometry.
  ExtrinsicGeometryInterface(HalfedgeMesh& mesh_);
  virtual ~ExtrinsicGeometryInterface() {}

public:
  /*
  // == Lengths, areas, and angles

  // Edge lengths
  inline void requireEdgeLengths();
  EdgeData<double> edgeLengths;

  // Face areas
  inline void requireFaceAreas();
  FaceData<double> faceAreas;
  */

protected:
};

} // namespace surface
} // namespace geometrycentral
