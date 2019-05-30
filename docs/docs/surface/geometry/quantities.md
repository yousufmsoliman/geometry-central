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

These quantities are defined for any `IntrinsicGeometry`, which is the base class of all other geometry objects---they will always be available on any kind of geometry.

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

    The "cotangent weight" of an interior halfedge, defined as $\frac{1}{2} \cot(\theta)$, where $\theta$ is the corner angle opposite the halfedge. Defined to be $0$ for exterior halfedges.

    Can be computed directly from edge lengths, or more efficiently in an embedded triangle via $\cot(\theta) = \frac{u \cdot v}{||u \times v||}$, where $u$ and $v$ are the edge vectors emanating from the opposite corner.

    Only valid on triangular meshes.

    - **member:** `HalfedgeData<double> IntrinsicGeometry::halfedgeCotanWeights`
    - **require:** `void IntrinsicGeometry::requireHalfedgeCotanWeights()`

??? func "edge cotan weight"

    The "cotangent weight" of an edge, defined as the sum of halfedge cotan weights for incident interior halfedges.

    Only valid on triangular meshes.

    - **member:** `EdgeData<double> IntrinsicGeometry::edgeCotanWeights`
    - **require:** `void IntrinsicGeometry::requireEdgeCotanWeights()`


## Tangent vectors and transport

These quantities are defined for any `IntrinsicGeometry`, which is the base class of all other geometry objects---they will always be available on any kind of geometry. Tangent vectors and transport are defined in terms of tangent spaces at faces and vertices, as defined below.

Recall that our `Vector2` types obey the multiplication and division rules of complex arithmetic, and thus can be used to represent rotations. For instance, a 2D vector representing a rotation can be used to rotate another vector like:
```cpp
Vector2 v = /* your vector */
Vector2 r = Vector2{std::cos(PI/4), std::sin(PI/4)}; // rotation by 45 degrees
Vector2 vRot = r * v;
```
This is fundamentally no different from using 2x2 rotation matrices, but leads to much cleaner code.

#### Face tangent spaces

To represent vectors that sit in flat mesh faces, we define a 2D coordinate frame tangent to each face. By default, this frame is aligned such that `face.halfedge()` points along the $x$-axis. All vectors in faces are then expressed via $(x,y)$ `Vector2D` coordinates in this frame. Crucially, this basis is well-defined even if the geometry does not have vertex positions.

#### Vertex tangent spaces

To represent vectors that sit at mesh faces, we consider a polar coordinate frame at each vertex. This frame is defined by measuring angles according to the rescaled corner angles as in `cornerScaledAngles`. By default, this frame is aligned such that `vertex.halfedge()` points along the $\phi=0$ $x$-axis. Of course, rather than using polar coordinates we can equivalently work in Cartesian frame---tangent vectors at vertices are then expressed via $(x,y)$ `Vector2D` coordinates in this frame. Crucially, this basis does not require picking a vertex normal, and is well-defined even if the geometry does not have vertex positions.

![vertex tangent coordinates diagram](../../../media/vertex_tangent_coordinates.svg)

??? func "halfedge vectors in face"

    Vectors for each halfedge in the coordinate frame of the face in which they sit. See the description of face tangent spaces above for a definition.

    Only valid on triangular meshes.

    - **member:** `HalfedgeData<Vector2> IntrinsicGeometry::halfedgeVectorsInFace`
    - **require:** `void IntrinsicGeometry::requireHalfedgeVectorsInFace()`


??? func "halfedge vectors in vertex"

    Vectors for each halfedge in the coordinate frame of the vertex from which the emanate (in `halfedge.vertex()`). See the description of vertex tangent spaces above for a definition.

    Only valid on triangular meshes.

    - **member:** `HalfedgeData<Vector2> IntrinsicGeometry::halfedgeVectorsInVertex`
    - **require:** `void IntrinsicGeometry::requireHalfedgeVectorsInVertex()`

??? func "transport vector across halfedge"

    Rotations which transport tangent vectors **across** a halfedge, rotating a vector from the tangent space of `halfedge.face()` to the tangent space `halfedge.twin().face()`.

    Always a unit vector, which can be multiplied by any other vector to compute the rotation. (recall our `Vector2`s multiply like complex numbers)

    Only valid on triangular meshes. Not defined for halfedges (interior or exterior) incident on boundary edges, these boundary values are set to NaN so errors can be caught quickly.

    - **member:** `HalfedgeData<Vector2> IntrinsicGeometry::transportVectorAcrossHalfedge`
    - **require:** `void IntrinsicGeometry::requireTransportVectorAcrossHalfedge()`
    
    Example usage:
    ```cpp
    geometry.requireTransportVectorAcrossHalfedge();

    Face f = /* ... */;        // a face of interest
    Vector2 myVec = /* ... */; // tangent vector in face f
    
    for(Halfedge he : f.adjacentHalfedges()) {

      Vertex neighborFace = he.twin().face();
      Vector2 rot = geometry.transportVectorAcrossHalfedge[he];
      Vector2 neighVec = rot * myVec;    // now in the basis of neighborFace
    }

    ```

??? func "transport vector along halfedge"

    Rotations which transport tangent vectors **along** a halfedge, rotating a vector from the tangent space of `halfedge.vertex()` to the tangent space `halfedge.twin().vertex()`.

    Always a unit vector, which can be multiplied by any other vector to compute the rotation. (recall our `Vector2`s multiply like complex numbers)

    Only valid on triangular meshes.

    - **member:** `HalfedgeData<Vector2> IntrinsicGeometry::transportVectorAlongHalfedge`
    - **require:** `void IntrinsicGeometry::requireTransportVectorAlongHalfedge()`
    
    Example usage:
    ```cpp
    geometry.requireTransportVectorAlongHalfedge();

    Vertex v = /* ... */;        // a vertex of interest
    Vector2 myVec = /* ... */;   // tangent vector in vertex v
    
    for(Halfedge he : v.outgoingHalfedges()) {
      Vertex neighborVertex = he.twin().vertex();
      Vector2 rot = geometry.transportVectorAlongHalfedge[he];
      Vector2 neighVec = rot * myVec;    // now in the basis of neighborVertex
    }

    ```


## Operators

## Extrinsic angles

## Embedded positions and normals
