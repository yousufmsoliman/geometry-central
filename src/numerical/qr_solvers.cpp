#include "geometrycentral/numerical/linear_solvers.h"

#include "geometrycentral/numerical/linear_algebra_utilities.h"

#ifdef HAVE_SUITESPARSE
#include "geometrycentral/numerical/suitesparse_utilities.h"
#include <SuiteSparseQR.hpp>
#include <cholmod.h>
#endif


using namespace Eigen;

namespace geometrycentral {

template <typename T>
Solver<T>::~Solver() {
#ifdef HAVE_SUITESPARSE
  if (cMat != nullptr) {
    cholmod_l_free_sparse(&cMat, context);
    cMat = nullptr;
  }
  if (cMatTrans != nullptr) {
    cholmod_l_free_sparse(&cMatTrans, context);
    cMatTrans = nullptr;
  }
  if (factorization != nullptr) {
    SuiteSparseQR_free(&factorization, context);
  }
#endif
}

template <typename T>
Solver<T>::Solver(SparseMatrix<T>& mat) : LinearSolver<T>(mat) {

// Check some sanity
#ifndef GC_NLINALG_DEBUG
  checkFinite(mat);
#endif

  mat.makeCompressed();

// Suitesparse version
#ifdef HAVE_SUITESPARSE

  // Is the system underdetermined?
  if (this->nRows < this-> : w gnCols) {
    underdetermined = true;
    throw std::logic_error("is not well tested, be careful");
  } else {
    underdetermined = false;
  }

  // Convert suitesparse format
  // Either use A or A^T, depending on whether the system underdetermined
  if (cMat != nullptr) {
    cholmod_l_free_sparse(&cMat, context);
  }
  cMat = toCholmod(this->mat, context);
  if (underdetermined) {
    if (cMatTrans != nullptr) {
      cholmod_l_free_sparse(&cMatTrans, context);
    }
    cMatTrans = cholmod_l_transpose(cMat, 2, context);
  }

  // Factor

  //  ordering options:
  //       0 or 3: fixed
  //       1: natural (only look for singletons)
  //       2: colamd after finding singletons
  //       4: CHOLMOD fter finding singletons
  //       5: amd(A'*A) after finding singletons
  //       6: metis(A'*A) after finding singletons
  //       7: SuiteSparseQR default (selects COLAMD, AMD, or METIS)
  const int ordering = 7;

  if (factorization != nullptr) {
    SuiteSparseQR_free(&factorization, context);
  }

  if (underdetermined) {
    factorization =
        SuiteSparseQR_factorize<typename Solver<T>::SOLVER_ENTRYTYPE>(ordering, zero_tolerance, cMatTrans, context);
  } else {
    factorization =
        SuiteSparseQR_factorize<typename Solver<T>::SOLVER_ENTRYTYPE>(ordering, zero_tolerance, cMat, context);
  }

  if (factorization == nullptr) {
    throw std::logic_error("Factorization failed");
  }

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
Vector<T> Solver<T>::solve(const Vector<T>& rhs) {
  Vector<T> out;
  solve(out, rhs);
  return out;
}

template <typename T>
void Solver<T>::solve(Vector<T>& x, const Vector<T>& rhs) {

  // Check some sanity
  if ((size_t)rhs.rows() != this->nRows) {
    throw std::logic_error("Vector is not the right length");
  }
#ifndef GC_NLINALG_DEBUG
  checkFinite(rhs);
#endif

// Suitesparse version
#ifdef HAVE_SUITESPARSE

  // Convert input to suitesparse format
  cholmod_dense* inVec = toCholmod(rhs, context);
  cholmod_dense* outVec;

  // Solve
  // outVec = SuiteSparseQR<typename Solver<T>::SOLVER_ENTRYTYPE>(cMat, inVec, context);
  // Note that the solve strategy is different for underdetermined systems
  if (underdetermined) {

    // solve y = R^-T b
    cholmod_dense* y =
        SuiteSparseQR_solve<typename Solver<T>::SOLVER_ENTRYTYPE>(SPQR_RTX_EQUALS_B, factorization, inVec, context);

    // compute x = Q*y
    outVec = SuiteSparseQR_qmult<typename Solver<T>::SOLVER_ENTRYTYPE>(SPQR_QX, factorization, y, context);
    cholmod_l_free_dense(&y, context);

  } else {

    // compute y = Q^T b
    cholmod_dense* y =
        SuiteSparseQR_qmult<typename Solver<T>::SOLVER_ENTRYTYPE>(SPQR_QTX, factorization, inVec, context);

    // solve x = R^-1 y
    // TODO what is this E doing here?
    outVec = SuiteSparseQR_solve<typename Solver<T>::SOLVER_ENTRYTYPE>(SPQR_RETX_EQUALS_B, factorization, y, context);

    cholmod_l_free_dense(&y, context);
  }

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
    // std::cerr << "Solver says: " << solver.lastErrorMessage() << std::endl;
    throw std::invalid_argument("Solve failed");
  }
#endif
}

template <typename T>
Vector<T> solve(SparseMatrix<T>& A, const Vector<T>& rhs) {
  Solver<T> s(A);
  return s.solve(rhs);
}


template <typename T>
size_t Solver<T>::rank() {
#ifdef HAVE_SUITESPARSE
  return factorization->rank;
#else
  return solver.rank();
#endif
}

// Explicit instantiations
template class Solver<double>;
template class Solver<float>;
template class Solver<std::complex<double>>;

template Vector<float> solve(SparseMatrix<float>& A, const Vector<float>& rhs);
template Vector<double> solve(SparseMatrix<double>& A, const Vector<double>& rhs);
template Vector<std::complex<double>> solve(SparseMatrix<std::complex<double>>& A,
                                            const Vector<std::complex<double>>& rhs);

} // namespace geometrycentral
