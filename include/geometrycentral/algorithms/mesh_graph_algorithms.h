#pragma once

#include "geometrycentral/geometry/geometry.h"
#include "geometrycentral/geometry/edge_length_geometry.h"
#include "geometrycentral/mesh/halfedge_data_types.h"

namespace geometrycentral {

// Find a subset of edges which connects all vertices
// Return value holds 'true' for an edge if it is in the tree
EdgeData<char> minimalSpanningTree(Geometry<Euclidean>* geometry);
EdgeData<char> minimalSpanningTree(EdgeLengthGeometry* geometry);

// Returns a set of edges which connect all vertices
// Note: Uses an MST+pruning approach to find short trees in O(N logN), but not guaranteed to be minimal; that's an
// NP-hard Steiner tree problem
EdgeData<char> spanningTreeBetweenVertices(Geometry<Euclidean>* geometry, const std::vector<VertexPtr>& requiredVertices);


} // namespace geometrycentral
