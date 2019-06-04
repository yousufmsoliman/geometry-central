#include "geometrycentral/surface/intrinsic_geometry_interface.h"

//#include "geometrycentral/surface/discrete_operators.h"

#include <fstream>
#include <limits>

namespace geometrycentral {
namespace surface {

// clang-format off
IntrinsicGeometryInterface::IntrinsicGeometryInterface(HalfedgeMesh& mesh_) : 
  BaseGeometryInterface(mesh_), 

  edgeLengthsQ              (&edgeLengths,          std::bind(&IntrinsicGeometryInterface::computeEdgeLengths, this),   quantities),
  faceAreasQ                (&faceAreas,            std::bind(&IntrinsicGeometryInterface::computeFaceAreas, this),     quantities)
  
  { }
// clang-format on

/*
void IntrinsicGeometryInterface::buildDependencies() {


  // === ALL the quantities
  // clang-format off
  //          quantity manager            dependencies                                  compute function
  //
  // == Basic geometric quantities
  addQuantity(faceAreasQ,                 {&edgeLengthsQ},                              &IntrinsicGeometryInterface::computeFaceAreas);
  addQuantity(vertexDualAreasQ,           {&faceAreasQ},                                &IntrinsicGeometryInterface::computeVertexDualAreas);
  addQuantity(edgeLengthsQ,               {},                                           &IntrinsicGeometryInterface::computeEdgeLengths);
  addQuantity(halfedgeCotanWeightsQ,      {&halfedgeOppositeAnglesQ},                   &IntrinsicGeometryInterface::computeHalfedgeCotanWeights);
  addQuantity(edgeCotanWeightsQ,          {&halfedgeCotanWeightsQ},                     &IntrinsicGeometryInterface::computeEdgeCotanWeights);
  addQuantity(vertexAngleDefectsQ,        {&halfedgeOppositeAnglesQ},                   &IntrinsicGeometryInterface::computeVertexAngleDefects);
  
  // == Vector fields, angles, and transport
  addQuantity(halfedgeFaceCoordsQ,              {&halfedgeOppositeAnglesQ, &edgeLengthsQ},          &IntrinsicGeometryInterface::computeHalfedgeFaceCoords);
  addQuantity(faceTransportCoefsQ,              {&halfedgeFaceCoordsQ},                             &IntrinsicGeometryInterface::computeFaceTransportCoefs);
  addQuantity(halfedgeOppositeAnglesQ,          {&edgeLengthsQ},                                    &IntrinsicGeometryInterface::computeHalfedgeOppositeAngles);
  addQuantity(halfedgeRescaledOppositeAnglesQ,  {&vertexAngleDefectsQ, &halfedgeOppositeAnglesQ},   &IntrinsicGeometryInterface::computeHalfedgeRescaledOppositeAngles);
  addQuantity(halfedgeVertexCoordsQ,            {&halfedgeRescaledOppositeAnglesQ},                 &IntrinsicGeometryInterface::computeHalfedgeVertexCoords);
  addQuantity(vertexTransportCoefsQ,            {&halfedgeVertexCoordsQ},                           &IntrinsicGeometryInterface::computeVertexTransportCoefs);
  
  // == Indices
  addQuantity(vertexIndicesQ,             {},                                           &IntrinsicGeometryInterface::computeVertexIndices);
  addQuantity(interiorVertexIndicesQ,     {},                                           &IntrinsicGeometryInterface::computeInteriorVertexIndices);
  addQuantity(faceIndicesQ,               {},                                           &IntrinsicGeometryInterface::computeFaceIndices);
  addQuantity(edgeIndicesQ,               {},                                           &IntrinsicGeometryInterface::computeEdgeIndices);
  addQuantity(halfedgeIndicesQ,           {},                                           &IntrinsicGeometryInterface::computeHalfedgeIndices);

  // == Operators
  addQuantity(basicDECOperatorsQ,         {&vertexDualAreasQ, &edgeCotanWeightsQ, &faceAreasQ, &vertexIndicesQ, &faceIndicesQ, &edgeIndicesQ},      &IntrinsicGeometryInterface::computeBasicDECOperators);
  addQuantity(zeroFormWeakLaplacianQ,     {&basicDECOperatorsQ},                        &IntrinsicGeometryInterface::computeZeroFormWeakLaplacian);
  // clang-format on
}
*/

// === Quantity implementations

void IntrinsicGeometryInterface::computeFaceAreas() {
  // ONEDAY try these for better accuracy in near-degenerate triangles?
  // "Miscalculating Area and Angles of a Needle-like Triangle" https://www.cs.unc.edu/~snoeyink/c/c205/Triangle.pdf

  faceAreas = FaceData<double>(mesh);
  for (Face f : mesh.faces()) {

    Halfedge he = f.halfedge();
    double a = edgeLengths[he.edge()];
    he = he.next();
    double b = edgeLengths[he.edge()];
    he = he.next();
    double c = edgeLengths[he.edge()];

    GC_SAFETY_ASSERT(he.next() == f.halfedge(), "faces mush be triangular"); 

    // Herons formula
    double s = (a + b + c) / 2.0;
    double area = std::sqrt(s * (s - a) * (s - b) * (s - c));

    faceAreas[f] = area;
  }
}


/*
void IntrinsicGeometryInterface::computeVertexDualAreas() {
  vertexDualAreas = VertexData<double>(*mesh);
  for (Vertex v : mesh->vertices()) {
    double A = 0;
    for (Face f : v.adjacentFaces()) {
      A += faceAreas[f];
    }
    vertexDualAreas[v] = A / 3.0;
  }
}

void IntrinsicGeometryInterface::computeHalfedgeFaceCoords() {
  verifyTriangular(mesh);
  halfedgeFaceCoords = HalfedgeData<Complex>(*mesh);

  for (Face f : mesh->faces()) {

    Halfedge he0 = f.halfedge();
    Halfedge he1 = he0.next();
    Halfedge he2 = he1.next();

    // Angles measured against he0
    halfedgeFaceCoords[he0] = Complex(edgeLengths[he0.edge()], 0.0);

    // Second halfedge
    double theta1 = PI - halfedgeOppositeAngles[he2];
    halfedgeFaceCoords[he1] = std::exp(IM_I * theta1) * Complex(edgeLengths[he1.edge()], 0.0);

    // Third halfedge
    double theta2 = halfedgeOppositeAngles[he1];
    halfedgeFaceCoords[he2] = -std::exp(IM_I * theta2) * Complex(edgeLengths[he2.edge()], 0.0);
  }

  for (Halfedge he : mesh->exteriorHalfedges()) {
    halfedgeFaceCoords[he] =
        std::numeric_limits<double>::quiet_NaN(); // using this basis is never a good idea, so NaN-out
  }
}


void IntrinsicGeometryInterface::computeFaceTransportCoefs() {

  faceTransportCoefs = HalfedgeData<Complex>(*mesh);

  for (Halfedge he : mesh->interiorHalfedges()) {
    if (he.twin().isInterior()) {
      Complex angleInSource = halfedgeFaceCoords[he];
      Complex desiredAngleInTarget = -halfedgeFaceCoords[he.twin()];
      faceTransportCoefs[he] = desiredAngleInTarget / angleInSource;
    }
  }
}

void IntrinsicGeometryInterface::computeVertexTransportCoefs() {

  vertexTransportCoefs = HalfedgeData<Complex>(*mesh);

  for (Halfedge he : mesh->halfedges()) {
    Complex angleInSource = halfedgeVertexCoords[he];
    Complex desiredAngleInTarget = -halfedgeVertexCoords[he.twin()];
    vertexTransportCoefs[he] = desiredAngleInTarget / angleInSource;
  }
}


void IntrinsicGeometryInterface::computeHalfedgeOppositeAngles() {
  verifyTriangular(mesh);

  halfedgeOppositeAngles = HalfedgeData<double>(*mesh);
  for (Halfedge he : mesh->halfedges()) {
    if (he.isInterior()) {

      double lOpp = edgeLengths[he.edge()];
      double lA = edgeLengths[he.next().edge()];
      double lB = edgeLengths[he.next().next().edge()];

      double q = (lA * lA + lB * lB - lOpp * lOpp) / (2. * lA * lB);
      q = clamp(q, -1.0, 1.0);
      double angle = std::acos(q);

      halfedgeOppositeAngles[he] = angle;
    } else {
      halfedgeOppositeAngles[he] = std::numeric_limits<double>::quiet_NaN();
    }
  }
}


void IntrinsicGeometryInterface::computeHalfedgeCotanWeights() {
  halfedgeCotanWeights = HalfedgeData<double>(*mesh);
  for (Halfedge he : mesh->halfedges()) {
    if (he.isInterior()) {
      halfedgeCotanWeights[he] = std::tan(PI / 2.0 - halfedgeOppositeAngles[he]);
    } else {
      halfedgeCotanWeights[he] = std::numeric_limits<double>::quiet_NaN();
    }
  }
}


void IntrinsicGeometryInterface::computeEdgeCotanWeights() {
  edgeCotanWeights = EdgeData<double>(*mesh);
  for (Edge e : mesh->edges()) {
    double weight = halfedgeCotanWeights[e.halfedge()];
    if (e.halfedge().twin().isInterior()) {
      weight += halfedgeCotanWeights[e.halfedge().twin()];
    }
    edgeCotanWeights[e] = 0.5 * weight;
  }
}


void IntrinsicGeometryInterface::computeVertexAngleDefects() {
  verifyTriangular(mesh);

  vertexAngleDefects = VertexData<double>(*mesh);
  for (Vertex v : mesh->vertices()) {
    if (v.isBoundary()) {
      // Convention: no curvature at boundary
      vertexAngleDefects[v] = 0;
    } else {
      double sum = 0;
      for (Halfedge he : v.outgoingHalfedges()) {
        sum += halfedgeOppositeAngles[he.next()];
      }
      vertexAngleDefects[v] = 2. * PI - sum;
    }
  }
}


void IntrinsicGeometryInterface::computeHalfedgeRescaledOppositeAngles() {
  halfedgeRescaledOppositeAngles = HalfedgeData<double>(*mesh);
  for (Halfedge he : mesh->interiorHalfedges()) {
    double origSum = 2. * PI - vertexAngleDefects[he.next().next().vertex()];
    halfedgeRescaledOppositeAngles[he] = halfedgeOppositeAngles[he] * 2. * PI / origSum;
  }
}


void IntrinsicGeometryInterface::computeHalfedgeVertexCoords() {
  verifyTriangular(mesh);

  halfedgeVertexCoords = HalfedgeData<Complex>(*mesh);

  for (Vertex v : mesh->vertices()) {

    if (v.isBoundary()) {

      // First, check what angle we associated with the boundary wedge
      // (recall that in a manifold triangle mesh, there can be at most one boundary wedge)
      double angleSum = 0;
      Halfedge afterBoundaryHe;
      for (Halfedge he : v.outgoingHalfedges()) {
        if (he.isInterior()) {
          angleSum += halfedgeRescaledOppositeAngles[he.next()];
        }
        if (!he.twin().isInterior()) {
          afterBoundaryHe = he;
        }
      }
      double boundaryAngle = 2 * PI - angleSum;

      // Now, loop like in the usual case, but substitute the boundary value when needed
      double coordSum = 0.0;

      // Custom loop to orbit CCW
      Halfedge firstHe = v.halfedge();
      Halfedge currHe = firstHe;
      do {
        halfedgeVertexCoords[currHe] = std::exp(coordSum * IM_I);
        if (currHe.isInterior()) {
          coordSum += halfedgeRescaledOppositeAngles[currHe.next()];
          currHe = currHe.next().next().twin();
        } else {
          coordSum += boundaryAngle;
          currHe = afterBoundaryHe;
        }
      } while (currHe != firstHe);


    } else {
      double coordSum = 0.0;

      // Custom loop to orbit CCW
      Halfedge firstHe = v.halfedge();
      Halfedge currHe = firstHe;
      do {
        halfedgeVertexCoords[currHe] = std::exp(coordSum * IM_I);
        coordSum += halfedgeRescaledOppositeAngles[currHe.next()];
        currHe = currHe.next().next().twin();
      } while (currHe != firstHe);
    }
  }
}


void IntrinsicGeometryInterface::computeBasicDECOperators() {

  { // Hodge 0
    size_t nVerts = mesh->nVertices();
    Eigen::VectorXd hodge0V(nVerts);
    for (Vertex v : mesh->vertices()) {
      double primalArea = 1.0;
      double dualArea = vertexDualAreas[v];
      double ratio = dualArea / primalArea;
      size_t iV = vertexIndices[v];
      hodge0V[iV] = ratio;
    }

    hodge0 = hodge0V.asDiagonal();
    hodge0Inv = hodge0V.asDiagonal().inverse();
  }


  { // Hodge 1
    size_t nEdges = mesh->nEdges();
    Eigen::VectorXd hodge1V(nEdges);
    for (Edge e : mesh->edges()) {
      double ratio = edgeCotanWeights[e];
      size_t iE = edgeIndices[e];
      hodge1V[iE] = ratio;
    }

    hodge1 = hodge1V.asDiagonal();
    hodge1Inv = hodge1V.asDiagonal().inverse();
  }

  { // Hodge 2
    size_t nFaces = mesh->nFaces();
    Eigen::VectorXd hodge2V(nFaces);
    for (Face f : mesh->faces()) {
      double primalArea = faceAreas[f];
      double dualArea = 1.0;
      double ratio = dualArea / primalArea;

      size_t iF = faceIndices[f];
      hodge2V[iF] = ratio;
    }
    hodge2 = hodge2V.asDiagonal();
    hodge2Inv = hodge2V.asDiagonal().inverse();
  }


  d0 = buildDerivative0(mesh);
  d1 = buildDerivative1(mesh);
}


void IntrinsicGeometryInterface::computeZeroFormWeakLaplacian() { zeroFormWeakLaplacian = d0.transpose() * hodge1 * d0; }
*/


} // namespace surface
} // namespace geometrycentral
