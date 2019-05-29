While the [halfedge mesh](../halfedge_mesh/basics.md) encodes the _connectivity_ of a surface, this section covers the classes which sit atop a halfedge mesh to define its _geometry_.

## Geometry hierarchy 

Geometry central is intentionally designed to allow flexibility in defining the geometry of a surface. Traditional algorithms typically define the geometry of a polygon mesh via a position for every vertex.

![geometry inheritance diagram](../../media/geometry_inheritance.svg)

## Cached quantities

We use a system of caches...
