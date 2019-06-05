#pragma once

#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"

#include <Eigen/SparseCore>

namespace geometrycentral {
namespace surface {


class EdgeLengthGeometry : public IntrinsicGeometryInterface {

public:
  EdgeLengthGeometry(HalfedgeMesh& mesh_, EdgeData<double>& inputEdgeLengths);
  virtual ~EdgeLengthGeometry() {}

  EdgeData<double> inputEdgeLengths;


protected:
  // Override the compute edge lengths method from intrinsic geometry.
  virtual void computeEdgeLengths() override;


private:
};

} // namespace surface
} // namespace geometrycentral
