#pragma once

#include "geometrycentral/geometry/dependent_quantity.h"
#include "geometrycentral/mesh/halfedge_mesh.h"
#include "geometrycentral/geometry/intrinsic_geometry.h"
#include "geometrycentral/utilities/unit_vector3.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

#include <Eigen/SparseCore>

#include <iostream>

namespace geometrycentral {


class EdgeLengthGeometry : public IntrinsicGeometry {

public:
  EdgeLengthGeometry(HalfedgeMesh* mesh_, EdgeData<double>& edgeLengths);
  EdgeLengthGeometry(HalfedgeMesh* mesh_, VertexData<Vector3>& vertexPositions);
  virtual ~EdgeLengthGeometry();

  void update(EdgeData<double> edgeLengths);

protected:
  std::vector<DependentQuantity*> allQuantities;

  // === Internal interface for all quantities

  DependentQuantity edgeLengthsQ;
  virtual void computeEdgeLengths() override;


private:
  EdgeData<double> geodesicEdgeLengths;
};

} // namespace geometrycentral
