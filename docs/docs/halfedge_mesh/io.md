## Reading meshes

Construct a halfedge mesh from a file on disk.

Include `geometrycentral/meshio.h`.

??? func "`#!cpp std::tuple<std::unique_ptr<HalfedgeMesh>,std::unique_ptr<Geometry>> loadMesh(std::string filename, std::string type="")`"

    Load a mesh from file. Returns both a `HalfedgeMesh` representing the connectivity, and a `Geometry` representing the geometry. See example below to concisely unpack.

    The `type` parameter determines the type of file to load. For example, `type="ply"` will attempt to read the target file as a .ply file. If no type is given, the type will be inferred from the file name. 

    Currently the following types are supported:
    
    - `obj`
    - `ply` (using [hapPLY](https://github.com/nmwsharp/happly))

    Example usage:
    ```cpp
    std::unique_ptr<HalfedgeMesh> mesh;
    std::unique_ptr<Geometry> geometry;
    std::tie<mesh, geometry> = loadMesh("spot.obj"); 
    ```

## Writing meshes

## Serializing containers 


Data stored in `MeshData<>` containers can be automatically written and loaded from file. Internally, data is stored as additional custom fields of a `.ply` file, though other software may not automatically understand these fields.

The `PlyHalfedgeMeshData` class is used to read and write these souped-up `.ply` files, and is distinct from the simple mesh-loading `.ply` interface above.

Include `geometrycentral/ply_halfedge_mesh_data.h`.

Example usage:


??? func "`#!cpp PlyHalfedgeMeshData::PlyHalfedgeMeshData(std::string filename, bool verbose=false)`"

    Construct from an existing `.ply` file.

??? func "`#!cpp PlyHalfedgeMeshData::PlyHalfedgeMeshData(HalfedgeMesh& mesh)`"

    Construct from an existing mesh. The mesh connectivity will be included when writing the file.

??? func "`#!cpp void PlyHalfedgeMeshData::write(std::string filename)`"

    Write the object to file.

### Writing properties

Add properties to the `PlyHalfedgeMeshData` object, which will be written when `write()` is called. The set of scalar types supported is the same as the [.ply](https://github.com/nmwsharp/happly) format, including list types.


??? func "`#!cpp void PlyHalfedgeMeshData::addVertexProperty<>(std::string name, VertexData<>& data)`"
    Add a property for writing. 

    - `name` A user-defined name with which the property will be written to file. Must be unique within the file.
    - `data` The data to be written, such a `VertexData<double>`.

??? func "`#!cpp void PlyHalfedgeMeshData::addHalfedgeProperty<>(std::string name, HalfedgeData<>& data)`"
    Add a property for writing. 

    - `name` A user-defined name with which the property will be written to file. Must be unique within the file.
    - `data` The data to be written, such a `HalfedgeData<double>`. 

??? func "`#!cpp void PlyHalfedgeMeshData::addEdgeProperty<>(std::string name, EdgeData<>& data)`"
    Add a property for writing. 

    - `name` A user-defined name with which the property will be written to file. Must be unique within the file.
    - `data` The data to be written, such a `EdgeData<double>`.

??? func "`#!cpp void PlyHalfedgeMeshData::addFaceProperty<>(std::string name, FaceData<>& data)`"
    Add a property for writing. 

    - `name` A user-defined name with which the property will be written to file. Must be unique within the file.
    - `data` The data to be written, such a `FaceData<double>`. 

### Reading properties
    
The template argument to this function will be likely necessary to resolve the expected type of the data. For instance, a property of type `double` on vertices could be accessed with.

```cpp
PlyHalfedgeMeshData data("my_file.ply");
VertexData<double> values = data.getVertexProperty<double>("propName");
```

Note that the automatic type promotion in [hapPLY](https://github.com/nmwsharp/happly) gives some flexibility in specifying the type. See the documentation there for details.

??? func "`#!cpp VertexData<T> PlyHalfedgeMeshData::getVertexProperty<T>(std::string name)`"
    Read a property from a loaded file. 

    - `name` A user-defined name with which the property will be read from the file. Throws if no such property exists.
    - **Return:** The requested container.

??? func "`#!cpp HalfedgeData<T> PlyHalfedgeMeshData::getHalfedgeProperty<T>(std::string name)`"
    Read a property from a loaded file. 

    - `name` A user-defined name with which the property will be read from the file. Throws if no such property exists.
    - **Return:** The requested container.

??? func "`#!cpp EdgeData<T> PlyHalfedgeMeshData::getEdgeProperty<T>(std::string name)`"
    Read a property from a loaded file. 

    - `name` A user-defined name with which the property will be read from the file. Throws if no such property exists.
    - **Return:** The requested container.

??? func "`#!cpp FaceData<T> PlyHalfedgeMeshData::getFaceProperty<T>(std::string name)`"
    Read a property from a loaded file. 

    - `name` A user-defined name with which the property will be read from the file. Throws if no such property exists.
    - **Return:** The requested container.



## Storing Delta-complexes

Most mesh file formats store connectivity via a face-vertex list; this format used automatically by the IO functions above. However, this format is insufficient for representing more general $\Delta$-complexes. To support IO for $\Delta$-complexes, connectivity can instead be encoded via halfedge adjacency indices as described in the [Internals](internals.md) section. This representation has the additional advantage that loading halfedge meshes will be very fast, as no connectivity needs to be detected.

The `.ply` readers automatically support reading this format. The option below enables writing `.ply` files in this format via the `PlyHalfedgeMeshData` class.


??? func "`#!cpp bool PlyHalfedgeMeshData::useHalfedgeAdjacency`"
    If true, writing will produce a `.ply` file which stores connectivity using haflfedge permutation indices rather than the usual face-vertex list.

    Default value: `false`.
