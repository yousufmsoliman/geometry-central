Generally, geometry central uses [Eigen](http://eigen.tuxfamily.org) for all matrix types. Though we build additional solvers and utilities on top of Eigen.

See [the Eigen section of dependencies](../../build/dependencies/#eigen) for instructions about getting Eigen and integrating with existing build systems.

Two typedefs are used extensively throughout geometry central to make Eigen types slightly less verbose. Both are defined in `linear_algebra_utilities.h`.

`#!cpp #include "geometrycentral/utilities/linear_algebra_utilities.h"`

??? func "`#!cpp Vector<T>`"

    A templated vector typedef, to Eigen's vector type.
    ```cpp
    template <typename T>
    using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    ```

    Use like `Vector<double>` or `Vector<bool>`.
    

??? func "`#!cpp SparseMatrix<T>`"

    A templated sparse matrix typedef, to Eigen's sparse matrix type.
    ```cpp
    template <typename T>
    using SparseMatrix = Eigen::SparseMatrix<T>;
    ```

    Use like `SparseMatrix<double>` or `SparseMatrix<int>`.
