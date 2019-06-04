#include "geometrycentral/surface/edge_length_geometry.h"

#include "geometrycentral/surface/discrete_operators.h"

#include <fstream>
#include <limits>

namespace geometrycentral {
namespace surface {

EdgeLengthGeometry::EdgeLengthGeometry(HalfedgeMesh& mesh_, EdgeData<double>& inputEdgeLengths_)
    : IntrinsicGeometryInterface(mesh_), inputEdgeLengths(inputEdgeLengths_) {}

void EdgeLengthGeometry::computeEdgeLengths() { edgeLengths = inputEdgeLengths; }

} // namespace surface
} // namespace geometrycentral
