#pragma once

#include "geometrycentral/surface/geometry.h"

#include <vector>


// Compute useful geometric quantities relating to optimal (Levi-Civita) transport

// Note: All of these functions implicitly follow the convention that angles in tangent space are measured against
// vertex.halfedge.

namespace geometrycentral {
namespace surface {

// === Completely compute direction fields
//     If the mesh has boundary, imposes dirichlet boundary conditions to
//     conform to the boundary.
//     Otherwise, computes the unit-norm solution
//     t \in [0,1] controls the strength of alignment with principal directions

VertexData<Vector2> computeSmoothestVertexDirectionField(IntrinsicGeomeryInterface& geometry, int nSym = 1,
                                                         bool alignCurvature = false);

FaceData<Vector2> computeSmoothestFaceDirectionField(IntrinsicGeomeryInterface& geometry, int nSym = 1,
                                                     bool alignCurvature = false);

// Find singularities in direction fields
FaceData<int> computeFaceIndex(IntrinsicGeomeryInterface& geometry, const VertexData<Vector2>& directionField, int nSym = 1);
VertexData<int> computeVertexIndex(IntrinsicGeomeryInterface& geometry, const FaceData<Complex>& directionField, int nSym = 1);

} // namespace surface
} // namespace geometrycentral
