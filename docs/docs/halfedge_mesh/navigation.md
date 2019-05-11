### Collection Iterators 

Use these routines to iterate over all of the elements in the mesh.

!!! note 
    Conveniently, these routines will work as expected if elements are _added_ to the mesh during iteration.

    However, _removing_ elements from the mesh in the midst of iteration, or calling `compress()` is not supported.

---

??? func "`#!cpp HalfedgeMesh::vertices()`"
    Iterate over the vertices in a mesh.
    ```cpp
    for(VertexPtr v : mesh.vertices()) {
      // do science here
    }
    ```

??? func "`#!cpp HalfedgeMesh::halfedges()`"
    Iterate over all of the halfedges in a mesh (both real and imaginary, if the mesh has boundary).
    ```cpp
    for(HalfedgePtr he : mesh.halfedges()) {
      // do science here
    }
    ```

??? func "`#!cpp HalfedgeMesh::realHalfedges()`"
    Iterate over the real halfedges in a mesh.
    ```cpp
    for(HalfedgePtr he : mesh.realHalfedges()) {
      // do science here
    }
    ```
    Note that on a boundary edge between vertices `i <--> j`, this set will only include a halfedge from `i --> j`, but not from `j --> i` (or vice versa).

??? func "`#!cpp HalfedgeMesh::imaginaryHalfedges()`"
    Iterate over the imaginary halfedges in a mesh.
    ```cpp
    for(HalfedgePtr he : mesh.imaginaryHalfedges()) {
      // do science here
    }
    ```
    Note that on a boundary edge between vertices `i <--> j`, this set will only include a halfedge from `i --> j`, but not from `j --> i` (or vice versa).

??? func "`#!cpp HalfedgeMesh::edges()`"
    Iterate over the edges in a mesh.
    ```cpp
    for(EdgePtr e : mesh.edges()) {
      // do science here
    }
    ```

??? func "`#!cpp HalfedgeMesh::faces()`"
    Iterate over the faces in a mesh.
    ```cpp
    for(FacePtr f : mesh.faces()) {
      // do science here
    }
    ```


## Neighborhood Iterators 

Use these routines to iterate over the neighbors of a mesh element.

---

!!! note
    The iterators in this section may have unexpected behavior in the case of a $\Delta$-complex, when there are (e.g.) self-edges, or multiple edges between a pair of vertices. Of course, for ordinary triangle mesh simplicial complexes they will behave as expected. See the [Delta complex](delta_complex.md) section for more information.

### Around a vertex

??? func "`#!cpp VertexPtr::outgoingHalfedges()`"
    Iterate over the halfedges which point outward from a vertex.
    ```cpp
    for(HalfedgePtr he : vert.outgoingHalfedges()) {
      assert(he.vertex() == vert); // true
      // do science here
    }
    ```

??? func "`#!cpp VertexPtr::incomingHalfedges()`"
    Iterate over the halfedges which point inward at a vertex.
    ```cpp
    for(HalfedgePtr he : vert.incomingHalfedges()) {
      assert(he.twin().vertex() == vert); // true
      // do science here
    }
    ```

??? func "`#!cpp VertexPtr::adjacentVertices()`"

    Iterate over the vertices edge-connected to this vertex.
    ```cpp
    for(VertexPtr v : vert.adjacentVertices()) {
      // do science here
    }
    ```

??? func "`#!cpp VertexPtr::adjacentEdges()`"

    Iterate over the edges incident on this vertex.
    ```cpp
    for(EdgePtr e : vert.adjacentEdges()) {
      // do science here
    }
    ```

??? func "`#!cpp VertexPtr::adjacentFaces()`"

    Iterate over the faces incident on this vertex.
    ```cpp
    for(FacePtr f : vert.adjacentFaces()) {
      // do science here
    }
    ```

### Around an edge

??? func "`#!cpp EdgePtr::adjacentHalfedges()`"

    Iterate over the two halfedges incident on this edge.
    ```cpp
    for(HalfedgePtr he : edge.adjacentHalfedges()) {
      // do science here
    }
    ```

??? func "`#!cpp EdgePtr::adjacentFaces()`"

    Iterate over the (one or two) faces incident on this edge.
    ```cpp
    for(FacePtr f : edge.adjacentFaces()) {
      // do science here
    }
    ```

### Around a face

??? func "`#!cpp FacePtr::adjacentVertices()`"
    Iterate over the vertices adjacent to a face.
    ```cpp
    for(VertexPtr v : face.adjacentVertices()) {
      // do science here
    }
    ```

??? func "`#!cpp FacePtr::adjacentHalfedges()`"
    Iterate over the halfedges incident on a face.
    ```cpp
    for(HalfedgePtr he : face.adjacentHalfedges()) {
      // do science here
    }
    ```

??? func "`#!cpp FacePtr::adjacentEdges()`"
    Iterate over the edges on the boundary of a face.
    ```cpp
    for(EdgePtr e : face.adjacentEdges()) {
      // do science here
    }
    ```

??? func "`#!cpp FacePtr::adjacentFaces()`"
    Iterate over the faces adjacent to a face, across each edge.
    ```cpp
    for(FacePtr f : face.adjacentFaces()) {
      // do science here
    }
    ```

## Accessors 

Use these routines to access elements of the mesh by their index.

!!! warning
    The indexing routines in the section are only valid when the mesh is _compressed_.

---

??? func "`#!cpp HalfedgePtr HalfedgeMesh::halfedge(size_t index)`"
    Constructs a reference to the i'th halfedge in the mesh. `0 <= index < nHalfedges()`.
    
    Warning: only valid when the mesh is _compressed_.

??? func "`#!cpp VertexPtr HalfedgeMesh::vertex(size_t index)`"
    Constructs a reference to the i'th vertex in the mesh. `0 <= index < nVertices()`.
    
    Warning: only valid when the mesh is _compressed_.

??? func "`#!cpp FacePtr HalfedgeMesh::face(size_t index)`"
    Constructs a reference to the i'th face in the mesh. `0 <= index < nFaces()`.
    
    Warning: only valid when the mesh is _compressed_.

??? func "`#!cpp EdgePtr HalfedgeMesh::edge(size_t index)`"
    Constructs a reference to the i'th edge in the mesh. `0 <= index < nEdges()`.
    
    Warning: only valid when the mesh is _compressed_.

??? func "`#!cpp FacePtr HalfedgeMesh::face(size_t index)`"
    Constructs a reference to the i'th face in the mesh. `0 <= index < nFaces()`.
    
    Warning: only valid when the mesh is _compressed_.
