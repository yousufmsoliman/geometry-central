#include "geometrycentral/surface/extrinsic_geometry_interface.h"

#include <limits>

namespace geometrycentral {
namespace surface {

// clang-format off
ExtrinsicGeometryInterface::ExtrinsicGeometryInterface(HalfedgeMesh& mesh_) : 
  IntrinsicGeometryInterface(mesh_)

  //edgeLengthsQ              (&edgeLengths,              std::bind(&ExtrinsicGeometryInterface::computeEdgeLengths, this),        {}),
  //faceAreasQ                (&faceAreas,                std::bind(&ExtrinsicGeometryInterface::computeFaceAreas, this),          {})
  
  {
    //quantities.push_back(&edgeLengthsQ);
    //quantities.push_back(&faceAreasQ);
  }
// clang-format on

} // namespace surface
} // namespace geometrycentral
