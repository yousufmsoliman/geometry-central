These routines allow modification of the mesh connectivity and insertion/deletion of elements.

As much as possible, these routines will check for validity before executing and throw an exception if something isn't right. The `NGC_SAFETY_CHECKS` define disables this behavior for a modest increase in performance, but is disabled by default even in release builds.

Note that aggressive use of these routines may reduce a mesh from a _simplicial complex_ to a _$\Delta$-complex_. For instance flipping enough edges in a mesh might create self-edges, which combinatorially connect a vertex to itself. See the [$\Delta$-complex](delta_complex.md) section for details, and an explanation of why these complexes are important.

<!--## Canonicalize-->
<!--**TODO** what canonical ordering even makes sense for the permutation-based mesh?-->

All operations invalidate the canonical ordering.

## Dynamic pointer types

A few of the operations listed below invalidate outstanding element references (like `HalfedgePtr`) by re-indexing the elements of the mesh. [Containers](containers.md) automatically update after re-indexing, and often code can be structured such that no element references need to be maintained across an invalidation.

However, if it is necessary to keep a reference to an element through a re-indexing, the `DynamicHalfedgePtr` can be used. These types behave like a `HalfedgePtr`, with the exception that they automatically update to remain valid when a mesh is re-indexed. These types should only be used when necessary, because they are expensive to maintain.

## In-place modifications

??? func "`#!cpp bool HalfedgeMesh::flip(EdgePtr e)`"

    Flip an edge by rotating counter-clockwise. 

    An edge cannot be combinatorially flipped if it is:

      - a boundary edge
      - incident on a degree-1 vertex.

    **Return:** true if the edge was actually flipped 

??? func "`#!cpp void HalfedgeMesh::setEdgeHalfedge(EdgePtr e, HalfedgePtr he)`"

    Re-index to ensure that `e.halfedge() == he`.
    
    `he` must be incident on `e`.

## Insertions

These routines modify a mesh by inserting new elements. Unlike deletions, insertions are totally transparent, and do not have any special consequences for other data structures and references. Element references remain valid, and [containers](containers.md) will automatically resize themselves to accommodate the new elements. 

Note that some operations my re-use existing elements to create their output. For instance, `splitEdge()` turns a single edge in to two; the index of the input edge will be re-used as one of the two output edges, and data along that edge will be unchanged in any containers.

!!! warning "Boundary loop invalidation"

    There is one special exception to the rule that inserting does not invalidate indexing. `FacePtr`s which point to boundary loops are invalidated after any operation which adds faces to the mesh. This is a consequence of the way we index boundary loops separate from faces, even though they are essentially faces in practice (see [Boundaries](boundaries.md) and [Internals](internals.md)) for details.

    As always, if a reference to a boundary loop must be preserved across an operation, a `DynamicFacePtr` will remain valid.

---

??? func "`#!cpp HalfedgePtr HalfedgeMesh::insertVertexAlongEdge(EdgePtr e)`"
    // Adds a vertex along an edge, increasing degree of faces. Returns ptr along the new edge, with he.vertex() as new
    // vertex and he.edge().halfedge() == he. Preserves canonical direction of edge.halfedge() for both halves of new
    // edge.
    HalfedgePtr insertVertexAlongEdge(EdgePtr e);

??? func "`#!cpp HalfedgePtr HalfedgeMesh::splitEdge(HalfedgePtr he)`"
    // Split an edge, also splitting adjacent faces. Returns new vertex.
    HalfedgePtr splitEdge(EdgePtr e);

??? func "`#!cpp HalfedgePtr HalfedgeMesh::splitEdge(EdgePtr e)`"

    Equivalent to `splitEdge(e.halfedge())`.

??? func "`#!cpp VertexPtr HalfedgeMesh::insertVertex(FacePtr f)`"

    // Add vertex inside face and triangulate. Returns new vertex.
    VertexPtr insertVertex(FacePtr f);


??? func "`#!cpp HalfedgePtr HalfedgeMesh::connectVertices(FacePtr face, VertexPtr vA, VertexPtr vB)`"

    // Same as above. Faster if you know the face.
    HalfedgePtr connectVertices(FacePtr face, VertexPtr vA, VertexPtr vB);

??? func "`#!cpp std::vector<FacePtr> HalfedgeMesh::triangulate(FacePtr face)`"

    // Triangulate in a face, returns all subfaces
    std::vector<FacePtr> triangulate(FacePtr face);


### Trimming storage
    
To amortize the cost of allocation, mesh buffers are resized sporadically in large increments; these resized buffers might significantly increase (e.g., double) the storage size of a mesh. Calling `trimStorage()` frees up any unused storage space to reduce memory usage. However, this function costs $\mathcal{O}(n)$ and should not be called in a tight loop.


??? func "`#!cpp void HalfedgeMesh::trimStorage()`"

    Free an additional storage associated with the mesh.

    As with insertions, does not invalidate references, except `FacePtr`s which point to boundary loops.


## Deletions

These routines delete mesh elements.


??? func "`#!cpp VertexPtr HalfedgeMesh::collapseEdge(EdgePtr e)`"

    // Collapse an edge. Returns the vertex adjacent to that edge which still exists. Returns VertexPtr() if not
    // collapsible.
    VertexPtr collapseEdge(EdgePtr e);

??? func "`#!cpp bool HalfedgeMesh::removeFaceAlongBoundary(FacePtr f)`"

    // Remove a face which is adjacent to the boundary of the mesh (along with its edge on the boundary).
    // Face must have exactly one boundary edge.
    // Returns true if could remove
    bool removeFaceAlongBoundary(FacePtr f);


### Compressed mode

Internally, the halfedge mesh is represented by dense arrays of indices which are lazily expanded (see [interals](internals.md) for details). To support fast deletion operations, we simply mark elements as deleted, without re-packing the index space. We say that the mesh is _compressed_ if the index space is dense and there are no such marked elements. When a mesh is not compressed, the `index` of a mesh element no longer serves as a proper enumeration from [0,N), but merely as a unique ID.

There are two consequences to being non-compressed:

  - Some operations cannot be implemented efficiently/correctly (e.g., random access of the i'th vertex)
  - Storage space is wasted by deleted elements


**All meshes are compressed after construction, and only become non-compressed if the user performs a deletion operation.**  The `makeCompressed()` function can be called to re-index the elements of the mesh as a proper enumeration from [0,N).

The `makeCompressed()` function invalidates pointers, and incurs an update of existing containers. As such, it is recommended to be called sporadically, after a sequence of operations is completed.

??? func "`#!cpp bool HalfedgeMesh::isCompressed()`"

    Returns true if the mesh is compressed.

??? func "`#!cpp void HalfedgeMesh::makeCompressed()`"

    Re-index the elements of the mesh to yield a dense enumeration. Invalidates all VertexPtr (etc) objects.

    Does nothing if the mesh is already compressed.

