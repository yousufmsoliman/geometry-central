These routines allow modifcation of the mesh connectivity.

As much as possible, these routines will check for validity before executing and throw an exception if something isn't right. The `NGC_SAFETY_CHECKS` define disables this behavior for a modest increase in performance, but is disabled by default even in release builds.

Note that agressive use of these routines may reduce a mesh from a _simplicial complex_ to a _$\Delta$-complex_. For instance flipping edges in a mesh might create self-edges, which combinatorially connect a vertex to itself. See the [$\Delta$-complex](delta_complex.md) section for details, and an explanation of why these complexes are important.

## Canonicalize

**TODO** what canonical ordering even makes sense for the permutation-based mesh?

All operations invalidate the canonical ordering.

## In-place modifications

??? func "`#!cpp bool HalfedgeMesh::flip(EdgePtr e)`"

    Flip an edge by rotating counter-clockwise. 

    An edge cannot be combinatorially flipped if it is:

      - a boundary edge
      - incident on a degree-1 vertex.

    **Return:** true if the edge was actually flipped 

??? func "`#!cpp void HalfedgeMesh::setEdgeHalfedge(EdgePtr e, HalfedgePtr he)`"

    Ensure that `e.halfedge() == he`.
    
    `he` must be incident on `e`.

## Insertions

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

??? func "`#!cpp void HalfedgeMesh::trimStorage()`"

    To amortize the cost of allocation, mesh buffers are resized sporadically in large increments; these resized buffers might significantly increase (e.g., double) the storage size of a mesh.

    Calling `trimStorage()` frees up any unused storage space, but ensures that the 


## Deletions
  
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

  - Some operations cannot be efficiently/correctly (e.g., looking up the i'th vertex)
  - Storage space is wasted by deleted elements

The `makeCompressed()` function can be called to re-index the elements of the mesh with a proper enumeration from [0,N).

Because the `makeCompressed()` function invalidates pointers an requires expensive 

??? func "`#!cpp void HalfedgeMesh::makeCompressed()`"

    Re-index the elements of the mesh to yield a dense enumeration. Invalidates all VertexPtr (etc) objects.

