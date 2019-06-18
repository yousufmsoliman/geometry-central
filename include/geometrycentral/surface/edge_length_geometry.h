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

  // == Immediates
  double faceArea(Face f) const;
  double cornerAngle(Corner c) const;
  double halfedgeCotanWeight(Halfedge he) const;
  double edgeCotanWeight(Edge e) const;


protected:
  // Override the compute edge lengths method from intrinsic geometry.
  virtual void computeEdgeLengths() override;


private:
};

} // namespace surface
} // namespace geometrycentral

#include "geometrycentral/surface/edge_length_geometry.ipp"
