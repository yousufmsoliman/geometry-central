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

  // Compute mean edge length and set shortTime
  double meanEdgeLength = 0.;
  for (Edge e : mesh.edges()) {
    meanEdgeLength += geom.edgeLengths[e];
  }
  meanEdgeLength /= mesh.nEdges();
  shortTime = tCoef * meanEdgeLength * meanEdgeLength;

  // We always want the mass matrix
  massMat = geom.vertexGalerkinMassMatrix;


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

VertexData<double> VectorHeatMethodSolver::extendScalar(const std::vector<std::tuple<SurfacePoint, double>>& sources) {
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
      double w = facePoint.faceCoords.y;
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



} // namespace surface
} // namespace geometrycentral
