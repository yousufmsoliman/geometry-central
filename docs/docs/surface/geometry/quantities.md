This page enumerates the surface geometry quantities available in geometry central.

Recall that these quantities are each associated with a [geometry interface](geometry.md#geometry-hierarchy) specifying what can be computed from the given input data. Instantiating a geometry from data, like a `VertexPositionGeometry` extends these interfaces and gives access to all of the quantities therein.  Quantities should usually be accessed via [managed caches](geometry.md#managed-quantities).

Typical usage:
```cpp
#include "geometrycentral/surface/geometry.h"
#include "geometrycentral/surface/meshio.h"
using namespace geometrycentral::surface;

// Load a mesh and geometry from file
std::unique_ptr<HalfedgeMesh> mesh;
std::unique_ptr<VertexPositionGeometry> positionGeometry;
std::tie<mesh, positionGeometry> = loadMesh("spot.obj");

// For the sake of the example, bind to an interface that offers
// only the quantities which we will actually use below.
IntrinsicGeometry& geometry = *positionGeometry;

// populate the quantity
geometry.requireFaceAreas();

for(Face f : mesh->faces()) {

  // access managed array holding quantity
  double area = geometry.faceAreas[f];

  // immediate computation: generally discouraged
  area = geometry.computeFaceArea(f);
}
```

## Lengths, areas, and intrinsic angles

Defined for any `IntrinsicGeometry`, which is the base class of all other geometry objects---these quantities will always be available on any kind of geometry.

??? func "edge length"

    The length of an edge in the mesh, as a non-negative real number.

    - **member:** `EdgeData<double> IntrinsicGeometry::edgeLengths`
    - **require:** `void IntrinsicGeometry::requireEdgeLengths()`
    - **immediate:** `double IntrinsicGeometry::computeEdgeLength(Edge e)`

??? func "face area"

    The area of a face, as a non-negative real number.

    May be computed from edge lengths via Heron's formula, or from embedded vertex positions with a cross product.

    Only valid on triangular meshes.

    - **member:** `FaceData<double> IntrinsicGeometry::faceAreas`
    - **require:** `void IntrinsicGeometry::requireFaceAreas()`
    - **immediate:** `double IntrinsicGeometry::computeFaceAreas(Face f)`

??? func "vertex dual area"

    An area associated with each vertex, as a non-negative real number.

    Only valid on triangular meshes.

    Defined to be $1/3$ the sum of all adjacent face areas. The sum of all vertex dual areas is equal to the usual surface area of the mesh.

    - **member:** `VertexData<double> IntrinsicGeometry::vertexDualAreas`
    - **require:** `void IntrinsicGeometry::requireVertexDualAreas()`

??? func "corner angles"

    The angle between incident edges at each corner of a mesh.

    Only valid on triangular meshes.

    - **member:** `CornerData<double> IntrinsicGeometry::cornerAngles`
    - **require:** `void IntrinsicGeometry::requireCornerAngles()`
    - **immediate:** `double IntrinsicGeometry::computeCornerAngle(Corner c)`

??? func "corner scaled angles"

    The angle between incident edges at each corner of a mesh, linearly rescaled such that the angles around every vertex sum to $2 \pi$. At boundary vertices, no scaling will be performed.

    Only valid on triangular meshes.

    - **member:** `CornerData<double> IntrinsicGeometry::cornerScaledAngles`
    - **require:** `void IntrinsicGeometry::requireCornerScaledAngles()`

??? func "vertex angle sum"

    The sum of corner angles around a vertex.

    Only valid on triangular meshes.

    - **member:** `VertexData<double> IntrinsicGeometry::vertexAngleSums`
    - **require:** `void IntrinsicGeometry::requireVertexAngleSums()`

??? func "vertex Gaussian curvature"

    The [_Gaussian curvature_](https://en.wikipedia.org/wiki/Gaussian_curvature) $K$ at a vertex, defined via the angle defect $K_v = 2 \pi - \sum \theta_i$, where $\sum \theta_i$ is the `vertexAngleSum` as above.

    Should be interpreted as an _integrated_ Gaussian curvature, giving the total curvature in the neighborhood of the vertex. On a closed surface, the [Gauss-Bonnet theorem](https://en.wikipedia.org/wiki/Gauss%E2%80%93Bonnet_theorem) tells us that the sum of these Gaussian curvatures will be a topological constant given by $\sum_v K_v = 2 \pi \chi$, where $\chi$ is the [Euler characteristic](../halfedge_mesh/basics.md#properties) of the surface. On surfaces with boundary, the geodesic curvature of the boundary factors in.

    Only valid on triangular meshes.

    - **member:** `VertexData<double> IntrinsicGeometry::vertexGaussianCurvatures`
    - **require:** `void IntrinsicGeometry::requireVertexGaussianCurvatures()`

??? func "face Gaussian curvature"

    The [_Gaussian curvature_](https://en.wikipedia.org/wiki/Gaussian_curvature) $K$ at a face, defined via the rescaled angle defect in the face $K_f = \pi - \sum \tilde{\theta}_i$, where $\tilde{\theta}_i$ are the _rescaled_ corner angles (as in `cornerScaledAngles`) incident on the face.

    Should be interpreted as an _integrated_ Gaussian curvature, giving the total curvature inside of the face. A corresponding curvature-per-unit-area can be computed by dividing by the area of the face.

    On a closed surface, the [Gauss-Bonnet theorem](https://en.wikipedia.org/wiki/Gauss%E2%80%93Bonnet_theorem) tells us that the sum of these Gaussian curvatures will be a topological constant given by $\sum_f K_f = 2 \pi \chi$, where $\chi$ is the [Euler characteristic](../halfedge_mesh/basics.md#properties) of the surface. On surfaces with boundary, the geodesic curvature of the boundary factors in.

    Only valid on triangular meshes.

    - **member:** `FaceData<double> IntrinsicGeometry::faceGaussianCurvatures`
    - **require:** `void IntrinsicGeometry::requireFaceGaussianCurvatures()`

??? func "halfedge cotan weight"

    The "cotangent weight" of an interior halfedge, defined as $\frac{1}{2} \cot(\theta)$, where $\theta$ is the corner angle opposite the halfedge. Defined to be $0$ for exterior faces.

    Can be computed directly from edge lengths, or more efficiently in an embedded triangle via $\cot(\theta) = \frac{u \cdot v}{||u \times v||}$, where $u$ and $v$ are the edge vectors emanating from the opposite corner.

    Only valid on triangular meshes.

    - **member:** `HalfedgeData<double> IntrinsicGeometry::halfedgeCotanWeights`
    - **require:** `void IntrinsicGeometry::requireHalfedgeCotanWeights()`

??? func "edge cotan weight"

    The "cotangent weight" of an edge, defined as the sum of halfedge cotan weights for incident interior halfedges.

    Only valid on triangular meshes.

    - **member:** `EdgeData<double> IntrinsicGeometry::edgeCotanWeights`
    - **require:** `void IntrinsicGeometry::requireEdgeCotanWeights()`

??? func "halfedge intrinsic vectors"

    Vectors for each halfedge in the coordinate frame of the face in which they sit. Computed by laying out the face isometrically in the plane via its edge lengths.

    This isometric layout can be used compute many other intrinsic geometry quantities, since it is well-defined even with purely intrinsic geometries.

    Only valid on triangular meshes.

    - **member:** `HalfedgeData<Vector2> IntrinsicGeometry::halfedgeIntrinsicVectors`
    - **require:** `void IntrinsicGeometry::requireHalfedgeIntrinsicVectors()`


## Tangent vectors and transport

## Operators

## Extrinsic angles

## Embedded positions and normals
