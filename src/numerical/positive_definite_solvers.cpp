#include "geometrycentral/numerical/linear_solvers.h"

#include "geometrycentral/numerical/linear_algebra_utilities.h"

#ifdef HAVE_SUITESPARSE
#include "geometrycentral/numerical/suitesparse_utilities.h"
#endif

using namespace Eigen;
using std::cout;
using std::endl;

namespace geometrycentral {

template <typename T>
PositiveDefiniteSolver<T>::~PositiveDefiniteSolver() {
#ifdef HAVE_SUITESPARSE
  if (cMat != nullptr) {
    cholmod_l_free_sparse(&cMat, context);
    cMat = nullptr;
  }
  if (factorization != nullptr) {
    cholmod_l_free_factor(&factorization, context);
  }
#endif
}

template <typename T>
PositiveDefiniteSolver<T>::PositiveDefiniteSolver(SparseMatrix<T>& mat) : LinearSolver<T>(mat) {


  // Check some sanity
  if (this->nRows != this->nCols) {
    throw std::logic_error("Matrix must be square");
  }
  size_t N = this->nRows;
#ifndef GC_NLINALG_DEBUG
  checkFinite(mat);
  checkHermitian(mat);
#endif

  mat.makeCompressed();

  // Suitesparse version
#ifdef HAVE_SUITESPARSE

  // Convert suitesparse format
  if (cMat != nullptr) {
    cholmod_l_free_sparse(&cMat, context);
  }
  cMat = toCholmod(mat, context, SType::SYMMETRIC);

  // Factor
  context.setSimplicial(); // must use simplicial for LDLt
  context.setLDL();        // ensure we get an LDLt factorization
  factorization = cholmod_l_analyze(cMat, context);
  cholmod_l_factorize(cMat, factorization, context);

  // Eigen version
#else
  solver.compute(mat);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Solver factorization error: " << solver.info() << std::endl;
    throw std::invalid_argument("Solver factorization failed");
  }
#endif
};

template <typename T>
Vector<T> PositiveDefiniteSolver<T>::solve(const Vector<T>& rhs) {
  Vector<T> out;
  solve(out, rhs);
  return out;
}

template <typename T>
void PositiveDefiniteSolver<T>::solve(Vector<T>& x, const Vector<T>& rhs) {

  size_t N = this->nRows;

  // Check some sanity
  if ((size_t)rhs.rows() != N) {
    throw std::logic_error("Vector is not the right length");
  }
#ifndef GC_NLINALG_DEBUG
  checkFinite(rhs);
#endif


  // Suitesparse version
#ifdef HAVE_SUITESPARSE

  // Convert input to suitesparse format
  cholmod_dense* inVec = toCholmod(rhs, context);

  // Solve
  cholmod_dense* outVec = cholmod_l_solve(CHOLMOD_A, factorization, inVec, context);

  // Convert back
  toEigen(outVec, context, x);

  // Free
  cholmod_l_free_dense(&outVec, context);
  cholmod_l_free_dense(&inVec, context);

  // Eigen version
#else
  // Solve
  x = solver.solve(rhs);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Solver error: " << solver.info() << std::endl;
    throw std::invalid_argument("Solve failed");
  }
#endif
}

template <typename T>
Vector<T> solvePositiveDefinite(SparseMatrix<T>& A, const Vector<T>& rhs) {
  PositiveDefiniteSolver<T> s(A);
  return s.solve(rhs);
}


// Explicit instantiations
template class PositiveDefiniteSolver<double>;
template class PositiveDefiniteSolver<float>;
template class PositiveDefiniteSolver<std::complex<double>>;

template Vector<float> solvePositiveDefinite<float>(SparseMatrix<float>& A, const Vector<float>& rhs);
template Vector<double> solvePositiveDefinite<double>(SparseMatrix<double>& A, const Vector<double>& rhs);
template Vector<std::complex<double>>
solvePositiveDefinite<std::complex<double>>(SparseMatrix<std::complex<double>>& A,
                                            const Vector<std::complex<double>>& rhs);


} // namespace geometrycentral
