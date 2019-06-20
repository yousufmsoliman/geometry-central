# Linear algebra utilities

### Construct and convert

??? func "`#!cpp SparseMatrix<T> identityMatrix(size_t N)`"

    Construct and `N x N` identity matrix of the requested type.


??? func "`#!cpp void shiftDiagonal(SparseMatrix<T>& m, T shiftAmount = 1e-4)`"

    Shift the diagonal of matrix, by adding `A + shiftDiagonal * identityMatrix()`.

??? func "`#!cpp SparseMatrix<double> complexToReal(const SparseMatrix<std::complex<double>>& m)`"

    Convert an `N x M` complex matrix to a `2N x 2M` real matrix, expanding each complex component in to a `2 x 2` to evaluate the complex product.


??? func "`#!cpp Vector<double> complexToReal(const Vector<std::complex<double>>& v)`"

    Convert an length `N` complex vector to a length `2N` real vector, expanding each complex component in to a consecutive real and imaginary component.


### Check matrix properties

### Block decomposition
