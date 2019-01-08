#include "geometrycentral/linear_solvers.h"

#include "geometrycentral/linear_algebra_utilities.h"


using namespace Eigen;

// Helper

namespace {

template <typename T>
double norm(Vector<T>& x, Eigen::SparseMatrix<T>& massMatrix) {
  return std::sqrt(std::abs((x.transpose() * massMatrix * x)[0]));
}
template double norm(Vector<double>& x, SparseMatrix<double>& massMatrix);
template double norm(Vector<float>& x, SparseMatrix<float>& massMatrix);
template double norm(Vector<geometrycentral::Complex>& x, SparseMatrix<geometrycentral::Complex>& massMatrix);

template <typename T>
void normalize(Vector<T>& x, Eigen::SparseMatrix<T>& massMatrix) {
  double scale = norm(x, massMatrix);
  x /= scale;
}
template void normalize(Vector<double>& x, SparseMatrix<double>& massMatrix);
template void normalize(Vector<float>& x, SparseMatrix<float>& massMatrix);
template void normalize(Vector<geometrycentral::Complex>& x, SparseMatrix<geometrycentral::Complex>& massMatrix);
}

namespace geometrycentral {

template <typename T>
Vector<T> smallestEigenvectorPositiveDefinite(SparseMatrix<T>& energyMatrix, SparseMatrix<T>& massMatrix,
                                              size_t nIterations) {

  // TODO could implement a faster variant in the suitesparse case; as-is this does a copy-convert each iteration

  size_t N = energyMatrix.rows();
  PositiveDefiniteSolver<T> solver(energyMatrix);

  Vector<T> u = Vector<T>::Random(N);
  Vector<T> x = u;
  for (size_t iIter = 0; iIter < nIterations; iIter++) {

    // Solve
    solver.solve(x, massMatrix * u);

    // Re-normalize
    normalize(x, massMatrix);

    // Update
    u = x;
  }

  return x;
}

template <typename T>
std::vector<Vector<T>> smallestKEigenvectorsPositiveDefinite(Eigen::SparseMatrix<T>& energyMatrix,
                                                             Eigen::SparseMatrix<T>& massMatrix, size_t kEigenvalues,
                                                             size_t nIterations) {

  std::vector<Vector<T>> res;

  size_t N = energyMatrix.rows();
  PositiveDefiniteSolver<T> solver(energyMatrix);

  auto projectOutPreviousVectors = [&](Vector<T>& x) {
    for (Vector<T>& v : res) {
      T proj = (((x.transpose() * massMatrix * v)[0]));
      x -= v * proj;
    }
  };

  for (size_t kEig = 0; kEig < kEigenvalues; kEig++) {

    Vector<T> u = Vector<T>::Random(N);
    projectOutPreviousVectors(u);
    Vector<T> x = u;
    for (size_t iIter = 0; iIter < nIterations; iIter++) {
      // Solve
      solver.solve(x, massMatrix * u);

      projectOutPreviousVectors(x);
      normalize(x, massMatrix);

      // Update
      u = x;
    }

    res.push_back(x);
  }

  return res;
}

template <typename T>
Vector<T> smallestEigenvectorSquare(SparseMatrix<T>& energyMatrix, SparseMatrix<T>& massMatrix, size_t nIterations) {

  // TODO could implement a faster variant in the suitesparse case; as-is this does a copy-convert each iteration

  size_t N = energyMatrix.rows();
  SquareSolver<T> solver(energyMatrix);

  Vector<T> u = Vector<T>::Random(N);
  Vector<T> x = u;
  for (size_t iIter = 0; iIter < nIterations; iIter++) {

    // Solve
    solver.solve(x, massMatrix * u);

    // Re-normalize
    normalize(x, massMatrix);

    // Update
    u = x;
  }

  return x;
}

template <typename T>
Vector<T> largestEigenvector(Eigen::SparseMatrix<T>& energyMatrix, Eigen::SparseMatrix<T>& massMatrix,
                             size_t nIterations) {

  size_t N = massMatrix.rows();
  PositiveDefiniteSolver<T> solver(massMatrix);

  Vector<T> u = Vector<T>::Random(N);
  Vector<T> x = u;
  for (size_t iIter = 0; iIter < nIterations; iIter++) {
    solver.solve(x, energyMatrix * u);
    normalize(x, massMatrix);
    u = x;
  }

  return x;
}

// Explicit instantiations
template Vector<double> smallestEigenvectorPositiveDefinite(SparseMatrix<double>& energyMatrix,
                                                            SparseMatrix<double>& massMatrix, size_t nIterations);
template Vector<float> smallestEigenvectorPositiveDefinite(SparseMatrix<float>& energyMatrix,
                                                           SparseMatrix<float>& massMatrix, size_t nIterations);
template Vector<Complex> smallestEigenvectorPositiveDefinite(SparseMatrix<Complex>& energyMatrix,
                                                             SparseMatrix<Complex>& massMatrix, size_t nIterations);

template std::vector<Vector<float>> smallestKEigenvectorsPositiveDefinite(Eigen::SparseMatrix<float>& energyMatrix,
                                                                          Eigen::SparseMatrix<float>& massMatrix,
                                                                          size_t kEigenvalues, size_t nIterations);
template std::vector<Vector<double>> smallestKEigenvectorsPositiveDefinite(Eigen::SparseMatrix<double>& energyMatrix,
                                                                           Eigen::SparseMatrix<double>& massMatrix,
                                                                           size_t kEigenvalues, size_t nIterations);
template std::vector<Vector<Complex>> smallestKEigenvectorsPositiveDefinite(Eigen::SparseMatrix<Complex>& energyMatrix,
                                                                            Eigen::SparseMatrix<Complex>& massMatrix,
                                                                            size_t kEigenvalues, size_t nIterations);

template Vector<double> smallestEigenvectorSquare(SparseMatrix<double>& energyMatrix, SparseMatrix<double>& massMatrix,
                                                  size_t nIterations);
template Vector<float> smallestEigenvectorSquare(SparseMatrix<float>& energyMatrix, SparseMatrix<float>& massMatrix,
                                                 size_t nIterations);
template Vector<Complex> smallestEigenvectorSquare(SparseMatrix<Complex>& energyMatrix,
                                                   SparseMatrix<Complex>& massMatrix, size_t nIterations);

template Vector<double> largestEigenvector(SparseMatrix<double>& energyMatrix, SparseMatrix<double>& massMatrix,
                                           size_t nIterations);
template Vector<float> largestEigenvector(SparseMatrix<float>& energyMatrix, SparseMatrix<float>& massMatrix,
                                          size_t nIterations);
template Vector<Complex> largestEigenvector(SparseMatrix<Complex>& energyMatrix, SparseMatrix<Complex>& massMatrix,
                                            size_t nIterations);


} // namespace geometrycentral
