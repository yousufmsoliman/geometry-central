#include "geometrycentral/surface/vector_heat_method.h"


#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

namespace geometrycentral {
namespace surface {

VectorHeatMethodSolver::VectorHeatMethodSolver(IntrinsicGeometryInterface& geom_, double tCoef_)
    : tCoef(tCoef_), mesh(geom_.mesh), geom(geom_)

{
  geom.requireEdgeLengths();
  geom.requireVertexGalerkinMassMatrix();
  geom.requireVertexLumpedMassMatrix();

  // Compute mean edge length and set shortTime
  double meanEdgeLength = 0.;
  for (Edge e : mesh.edges()) {
    meanEdgeLength += geom.edgeLengths[e];
  }
  meanEdgeLength /= mesh.nEdges();
  shortTime = tCoef * meanEdgeLength * meanEdgeLength;

  // We always want the mass matrix
  // massMat = geom.vertexGalerkinMassMatrix; TODO
  massMat = geom.vertexLumpedMassMatrix;


  geom.unrequireVertexLumpedMassMatrix();
  geom.unrequireVertexGalerkinMassMatrix();
  geom.unrequireEdgeLengths();
}


void VectorHeatMethodSolver::ensureHaveScalarHeatSolver() {
  if (scalarHeatSolver != nullptr) return;

  // Get the ingredients
  geom.requireCotanLaplacian();
  SparseMatrix<double>& L = geom.cotanLaplacian;

  // Build the operator
  SparseMatrix<double> heatOp = massMat + shortTime * L;
  scalarHeatSolver.reset(new PositiveDefiniteSolver<double>(heatOp));

  geom.unrequireCotanLaplacian();
}

void VectorHeatMethodSolver::ensureHaveVectorHeatSolver() {
  if (vectorHeatSolver != nullptr) return;

  // Get the ingredients
  geom.requireVertexConnectionLaplacian();
  SparseMatrix<std::complex<double>>& Lconn = geom.vertexConnectionLaplacian;

  // Build the operator
  SparseMatrix<std::complex<double>> vectorOp = massMat.cast<std::complex<double>>() + shortTime * Lconn;
  vectorHeatSolver.reset(new PositiveDefiniteSolver<std::complex<double>>(vectorOp));

  geom.unrequireVertexConnectionLaplacian();
}


void VectorHeatMethodSolver::ensureHavePoissonSolver() {
  if (poissonSolver != nullptr) return;

  // Get the ingredients
  geom.requireCotanLaplacian();
  SparseMatrix<double>& L = geom.cotanLaplacian;

  // Build the operator
  poissonSolver.reset(new PositiveDefiniteSolver<double>(L));

  geom.unrequireCotanLaplacian();
}

VertexData<double> VectorHeatMethodSolver::extendScalar(const std::vector<std::tuple<Vertex, double>>& sources) {

  std::vector<std::tuple<SurfacePoint, double>> sourcePoints;
  for (auto tup : sources) {
    sourcePoints.emplace_back(SurfacePoint(std::get<0>(tup)), std::get<1>(tup));
  }

  // call general version
  return extendScalar(sourcePoints);
}

VertexData<double> VectorHeatMethodSolver::extendScalar(const std::vector<std::tuple<SurfacePoint, double>>& sources) {
  if (sources.size() == 0) {
    return VertexData<double>(mesh, std::numeric_limits<double>::quiet_NaN());
  }

  ensureHaveScalarHeatSolver();

  geom.requireVertexIndices();

  // === Build the RHS
  Vector<double> dataRHS = Vector<double>::Zero(mesh.nVertices());
  Vector<double> indicatorRHS = Vector<double>::Zero(mesh.nVertices());

  for (auto tup : sources) {
    SurfacePoint point = std::get<0>(tup);
    double value = std::get<1>(tup);

    SurfacePoint facePoint = point.inSomeFace();
    Halfedge he = facePoint.face.halfedge();

    { // First adjacent vertex
      size_t vInd = geom.vertexIndices[he.vertex()];
      double w = facePoint.faceCoords.x;
      dataRHS[vInd] += w * value;
      indicatorRHS[vInd] += w;
    }
    he = he.next();

    { // Second adjacent vertex
      size_t vInd = geom.vertexIndices[he.vertex()];
      double w = facePoint.faceCoords.y;
      dataRHS[vInd] += w * value;
      indicatorRHS[vInd] += w;
    }
    he = he.next();

    { // Third adjacent vertex
      size_t vInd = geom.vertexIndices[he.vertex()];
      double w = facePoint.faceCoords.z;
      dataRHS[vInd] += w * value;
      indicatorRHS[vInd] += w;
    }
  }


  // == Solve the systems
  Vector<double> dataSol = scalarHeatSolver->solve(dataRHS);
  Vector<double> indicatorSol = scalarHeatSolver->solve(indicatorRHS);


  // == Combine results
  Vector<double> interpResult = dataSol.array() / indicatorSol.array();
  VertexData<double> result(mesh, interpResult);

  geom.unrequireVertexIndices();

  return result;
}


VertexData<Vector2>
VectorHeatMethodSolver::transportTangentVectors(const std::vector<std::tuple<SurfacePoint, Vector2>>& sources) {
  if (sources.size() == 0) {
    return VertexData<Vector2>(mesh, Vector2::undefined());
  }

  geom.requireVertexIndices();


  // === Setup work

  // Don't need to do magnitude solve with a single source
  bool singleVec = sources.size() == 1;

  // Make sure systems have been built and factored
  ensureHaveVectorHeatSolver();
  if (!singleVec) {
    ensureHaveScalarHeatSolver();
  }


  // === Build the RHS

  Vector<std::complex<double>> dirRHS = Vector<std::complex<double>>::Zero(mesh.nVertices());

  // Accumulate magnitude data for scalar problem
  std::vector<std::tuple<SurfacePoint, double>> magnitudeSources;

  for (auto tup : sources) {
    SurfacePoint point = std::get<0>(tup);
    Vector2 vec = std::get<1>(tup);
    std::complex<double> unitVec = Vector2::fromComplex(vec).normalize();

    // Add to the list of magnitudes for magnitude interpolation
    magnitudeSources.emplace_back(point, vec.norm());

    SurfacePoint facePoint = point.inSomeFace();
    Halfedge he = facePoint.face.halfedge();

    { // First adjacent vertex
      size_t vInd = geom.vertexIndices[he.vertex()];
      double w = facePoint.faceCoords.x;
      dirRHS[vInd] += w * unitVec;
    }
    he = he.next();


    { // Second adjacent vertex
      size_t vInd = geom.vertexIndices[he.vertex()];
      double w = facePoint.faceCoords.y;
      dirRHS[vInd] += w * unitVec;
    }
    he = he.next();


    { // Third adjacent vertex
      size_t vInd = geom.vertexIndices[he.vertex()];
      double w = facePoint.faceCoords.z;
      dirRHS[vInd] += w * unitVec;
    }
  }

  auto vectorQ = polyscope::getSurfaceMesh()->addVertexIntrinsicVectorQuantity("rhs vectors", dirRHS);


  // == Solve the system

  Vector<std::complex<double>> vecSolution = vectorHeatSolver->solve(dirRHS);


  // == Get the magnitude right

  VertexData<Vector2> result(mesh);
  if (singleVec) {
    // For one sources, can just normalize and project
    double targetNorm = std::get<1>(sources[0]).norm();
    std::cout << "target norm = " << targetNorm << std::endl;

    // for (size_t i = 0; i < (size_t)vecSolution.rows(); i++) {
    // std::cout << "i = " << i << " sol = " << vecSolution[i] << " sol.abs() = " << std::abs(vecSolution[i])
    //<< " sol norm = " << vecSolution[i] / std::abs(vecSolution[i])
    //<< " sol norm mag = " << std::abs(vecSolution[i] / std::abs(vecSolution[i])) << std::endl;
    //}

    vecSolution = (vecSolution.array() / vecSolution.array().abs()) * targetNorm;
    // for (size_t i = 0; i < (size_t)vecSolution.rows(); i++) {
    // std::cout << "i = " << i << " sol = " << vecSolution[i] << " sol.abs() = " << std::abs(vecSolution[i])
    //<< std::endl;
    //}

    // Copy to output vector
    for (Vertex v : mesh.vertices()) {
      result[v] = Vector2::fromComplex(vecSolution[geom.vertexIndices[v]]);
      std::cout << "v = " << v << " sol = " << result[v] << " sol.abs() = " << norm(result[v]) << std::endl;
    }
  } else {
    // For multiple sources, need to interpolate magnitudes

    // === Perform scalar interpolation
    VertexData<double> interpMags = extendScalar(magnitudeSources);

    // Scale and copy to result
    for (Vertex v : mesh.vertices()) {
      Vector2 dir = Vector2::fromComplex(vecSolution[geom.vertexIndices[v]]).normalize();
      result[v] = dir * interpMags[v];
    }
  }


  geom.unrequireVertexIndices();
  return result;
}


} // namespace surface
} // namespace geometrycentral
