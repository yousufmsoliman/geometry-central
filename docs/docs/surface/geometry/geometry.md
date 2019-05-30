While the [halfedge mesh](../halfedge_mesh/basics.md) encodes the _connectivity_ of a surface, this section covers the classes which sit atop a halfedge mesh to define its _geometry_.

## Geometry hierarchy 


!!! TLDR "TL;DR"

    Construct a `VertexPositionGeometry` object using vertex positions; it offers all the geometric routines you would expect, and can be passed to any method that demands geometry.

    Many algorithms can actually operate on weaker data than vertex positions. Read on to learn more.



Geometry central is intentionally designed to allow flexibility in defining the geometry of a surface. Traditional code might assume a 3D position for every vertex, but many algorithms actually need only the _intrinsic geometry_ of a surface, aka the edge lengths. More generally, specifying algorithms to only use the geometric data they really need allows us to seamlessly leverage powerful techniques.

We (sparingly) make use of polymorphism via inheritance in C++ to encode a hierarchy of geometric quantities that one might compute for a surface. 

- **Interfaces** define which quantities can be computed from the geometry; for instance, an `EmbeddedGeometryInterface` can compute face normals, and it can also compute face areas because it extends the more basic `IntrinsicGeometryInterface`. Interfaces are abstract, and cannot be instantiated by themselves.
- **Realizations** are concrete classes allow the user instantiate a geometry object from data; for instance, a `VertexPositionGeometry` can be constructed from vertex positions, and implements the `EmbeddedGeometryInterface` giving access to a wide range of intrinsic and extrinsic geometric quantities.

The following diagram outlines the interfaces and realizations currently available.

![geometry inheritance diagram](../../media/geometry_inheritance.svg)

## Cached quantities

We use a system of caches...
