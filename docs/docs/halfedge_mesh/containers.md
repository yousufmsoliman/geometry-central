The halfedge mesh class is equipped with a system of containers for associating data with mesh vertices, halfedges, edges, and faces. For instance, to represent a scalar value at vertices, or a vector value at faces, one can use 

```cpp
// on vertices
VertexData<double> myVertexScalar(mesh);
VertexPtr v = /* some vertex */;
myVertexScalar[v] = 42.;

// on faces
FaceData<Vector3> myFaceVector(mesh);
FacePtr f = /* some face */;
myFaceVector[f] = Vector3{1., 2., 3.};
```
and so on.

## Mesh data types

The mesh data types are all templated on a common base class: `MeshData<P,T>`, where `P` is a pointer type (such as `VertexPtr`) and `T` is a scalar type (such as `double`). The first template should usually be omitted in user code; the various element containers are all typedef'd with concise names as follows:

- `VertexData<T>` data at vertices
- `HalfedgeData<T>` data at (interior and exterior) halfedges 
- `EdgeData<T>` data at edges 
- `FaceData<T>` data at faces 
- `BoundaryLoopData<T>` data at boundary loops

Most functionality is identical between all of these classes, so the sections below are written in terms of the generic `MeshData<>` class.

## Construction


??? func "`#!cpp MeshData<P,T>::MeshData<P,T>(HalfedgeMesh& mesh)`"
    Construct a new container over a mesh. Elements will be default-initialized with `T()`.

??? func "`#!cpp MeshData<P,T>::MeshData<P,T>(HalfedgeMesh& mesh, T initVal)`"
    Construct a new container over a mesh. 
    
    Elements will be default-initialized with `initVal`, and any newly-created mesh elements will have their default values set to `initVal`.

Additionally, see the vector-based initializers in [vector interop](containers.md#vector-interop).

## Accessors

## Mutation

## Vector interop

??? func "`#!cpp MeshData<P,T>::MeshData<P,T>(HalfedgeMesh& mesh)`"
    Construct a new container over a mesh. Elements will be default-initialized with `T()`.

## Advanced features



### Oriented edge data 

<!--TODO reword...-->
Scalar values on edges often carry meaning with respect to some oriented direction along the edge--- common examples include differences between values at vertices, or more generally 1-forms in discrete differential geometry. In such settings, a scalar value is concisely stored along edges, but its sign should flip when accessed "along" the opposite direction.

`EdgeData<T>` containers offer a pair of special additional accessors for oriented data, which handle the sign flips automatically.

??? func "`#!cpp T EdgeData<T>::getOriented(HalfedgePtr he)`"

    Access edge-valued data with sign determined by canonical halfedge orientation.

    Returns `edgeData[he.edge()]` if `he == he.edge().halfedge()`, or `-edgeData[he.edge()]` otherwise.


??? func "`#!cpp void EdgeData<T>::setOriented(HalfedgePtr he, T val)`"
    
    Access edge-valued data with sign determined by canonical halfedge orientation.

    Sets `edgeData[he.edge()] = val` if `he == he.edge().halfedge()`, or `edgeData[he.edge()] = -val` otherwise.
