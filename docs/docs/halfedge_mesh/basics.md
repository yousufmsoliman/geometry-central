# Halfedge meshes

The halfedge mesh is a powerful and flexible data structure for representing oriented, manifold polygonal meshes, and is the core data structure in geometry-central.

The halfedge mesh has several key advantages over other data structures, most notably that all adjacent-neighborhood traversals can be implemented in constant time, without the use of any variably-sized neighbor lists. Furthermore, common mutation operations like edge splits and vertex insertions can be performed in constant time.

As the name suggests, the primary type in a halfedge mesh is a _halfedge_, in addition to the usual _vertex_, _edge_ and _face_ types. A halfedge is a directed edge incident on a face, as shown below. Two halfedges, oriented in opposite directions, make up each edge in the mesh. Each halfedge has relationships with five adjacent elements: 

- `HalfedgePtr::twin()` the other halfedge across the incident edge
- `HalfedgePtr::next()` the next halfedge in clockwise order around the incident face
- `HalfedgePtr::vertex()` the vertex at the tail (back) of the halfedge
- `HalfedgePtr::edge()` the incident edge
- `HalfedgePtr::face()` the incident face

![halfedge pointers](../media/halfedge_pointers.png)

Each vertex, edge, and face need just one relationship:

- `VertexPtr::halfedge()` _any_ of the incident halfedges (which point outward from the vertex)
- `EdgePtr::halfedge()` _any_ of the incident halfedges
- `FacePtr::halfedge()` _any_ of the incident halfedges

Notice how this fixed set of relationships can be used implement local traversals. For instance, the neighboring face across an edge can be accessed with `face.halfedge().twin().face()`, and we can iterate around the neighbors of a vertex by advancing `currHe = currHe.twin().next()` and examining `currHe.vertex()`. See [navigation](navigation.md) for more information on traversals, and convenience iterators.

You many notice in the above examples that the primary type we use to interact with halfedge mesh elements is the lightweight `HalfedgePtr` (etc) types. These types logically refer to a mesh element, and offer routines to access their neighbors. There is no explicit `Halfedge` class in our API, only references to logical halfedges.

## Manifold, Oriented Surfaces

The basic halfedge mesh imposes two requirements: manifoldness and orientability. 

Manifoldness means that our surface must locally look like a plane in any neighborhood. This disallows structures such as three faces meeting at an edge, or two cones of faces meeting at a single vertex like an hourglass. 

Furthermore the halfedge mesh implies a combinatorial _orientation_ of the surface, indicated by the clockwise ordering of halfedges around each face (see figure below). Because the halfedge mesh implies an orientation, it cannot represent non-orientable surfaces, like a Klein bottle.

![halfedge orientation](../media/halfedge_orientation.png)

These properties are invariants which always hold for any meaningful halfedge mesh; in practice we check them during construction and ensure that all operations preserve them.

Note that our halfedge mesh _does not_ require that faces be triangles or quads; arbitrary faces with degree >= 3 are supported. However, many operations are only defined for triangle meshes and will throw errors if invoked on other meshes.

## Boundaries



## Basic API


### Constructors

### Element counts

??? func "`#!cpp size_t HalfedgeMesh::nVertices()`"
    Returns the number of vertices. 

??? func "`#!cpp size_t HalfedgeMesh::nInteriorVertices()`"
    Returns the number of vertices not incident on the boundary.

??? func "`#!cpp size_t HalfedgeMesh::nBoundaryVertices()`"
    Returns the number of vertices incident on the boundary.

??? func "`#!cpp size_t HalfedgeMesh::nEdges()`"
    Returns the number of edges. 

??? func "`#!cpp size_t HalfedgeMesh::nFaces()`"
    Returns the number of faces. Does not include imaginary "faces" created around boundaries.

??? func "`#!cpp size_t HalfedgeMesh::nHalfedges()`"
    Returns the number of halfedges, including both real and halfedges and any imaginary halfedges incident on boundary loops. Always exactly twice the number of edges.

??? func "`#!cpp size_t HalfedgeMesh::nRealHalfedges()`"
    Returns the number of real halfedges, which are incident on faces of the mesh. Always equal to the sum of the number of sides of all faces.

??? func "`#!cpp size_t HalfedgeMesh::nImaginaryHalfedges()`"
    Returns the number of imaginary halfedges, which are opposite boundary faces. 

??? func "`#!cpp size_t HalfedgeMesh::nBoundaryLoops()`"
    Returns the number of distinct boundary loops in the mesh, each identified as an "imaginary" face closing a boundary loop in the mesh.


### Properties

??? func "`#!cpp int HalfedgeMesh::eulerCharacteristic()`"
    Returns the Euler characteristic of the surface. Computed in O(1) from element counts. Warning: assumes the mesh is a single connected component.

??? func "`#!cpp int HalfedgeMesh::genus()`"
    Returns the genus of the surface. Computed in O(1) from element counts. Warning: assumes the mesh is a single connected component.

??? func "`#!cpp bool HalfedgeMesh::isTriangular()`"
    Returns true if all faces in the mesh have 3 sides. Complexity O(n), do not call in a tight loop.

??? func "`#!cpp size_t HalfedgeMesh::nConnectedComponents()`"
    Returns the number of distinct connected components of the . Complexity O(n), do not call in a tight loop.


### Utility functions   

??? func "`#!cpp std::unique_ptr<Halfedgemesh> copy()`"
    Constructs a copy of the mesh.


