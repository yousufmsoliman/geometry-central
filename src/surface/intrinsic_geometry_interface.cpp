#include "geometrycentral/surface/intrinsic_geometry_interface.h"

//#include "geometrycentral/surface/discrete_operators.h"

#include <fstream>
#include <limits>

using std::cout;
using std::endl;

namespace geometrycentral {
namespace surface {

// clang-format off
IntrinsicGeometryInterface::IntrinsicGeometryInterface(HalfedgeMesh& mesh_) : 
  BaseGeometryInterface(mesh_), 

  edgeLengthsQ              (&edgeLengths,                  std::bind(&IntrinsicGeometryInterface::computeEdgeLengths, this),               quantities),
  faceAreasQ                (&faceAreas,                    std::bind(&IntrinsicGeometryInterface::computeFaceAreas, this),                 quantities),
  vertexDualAreasQ          (&vertexDualAreas,              std::bind(&IntrinsicGeometryInterface::computeVertexDualAreas, this),           quantities),
  cornerAnglesQ             (&cornerAngles,                 std::bind(&IntrinsicGeometryInterface::computeCornerAngles, this),              quantities),
  vertexAngleSumsQ          (&vertexAngleSums,              std::bind(&IntrinsicGeometryInterface::computeVertexAngleSums, this),           quantities),
  cornerScaledAnglesQ       (&cornerScaledAngles,           std::bind(&IntrinsicGeometryInterface::computeCornerScaledAngles, this),        quantities),
  vertexGaussianCurvaturesQ (&vertexGaussianCurvatures,     std::bind(&IntrinsicGeometryInterface::computeVertexGaussianCurvatures, this),  quantities),
  faceGaussianCurvaturesQ   (&faceGaussianCurvatures,       std::bind(&IntrinsicGeometryInterface::computeFaceGaussianCurvatures, this),    quantities),
  halfedgeCotanWeightsQ     (&halfedgeCotanWeights,         std::bind(&IntrinsicGeometryInterface::computeHalfedgeCotanWeights, this),      quantities),
  edgeCotanWeightsQ         (&edgeCotanWeights,             std::bind(&IntrinsicGeometryInterface::computeEdgeCotanWeights, this),          quantities)
  
  { }
// clang-format on

/*
void IntrinsicGeometryInterface::buildDependencies() {


  // === ALL the quantities
  // clang-format off
  //          quantity manager            dependencies                                  compute function
  //
  // == Basic geometric quantities
  addQuantity(faceAreasQ,                 {&edgeLengthsQ}, &IntrinsicGeometryInterface::computeFaceAreas);
  addQuantity(vertexDualAreasQ,           {&faceAreasQ}, &IntrinsicGeometryInterface::computeVertexDualAreas);
  addQuantity(edgeLengthsQ,               {}, &IntrinsicGeometryInterface::computeEdgeLengths);
  addQuantity(halfedgeCotanWeightsQ,      {&halfedgeOppositeAnglesQ},
&IntrinsicGeometryInterface::computeHalfedgeCotanWeights); addQuantity(edgeCotanWeightsQ, {&halfedgeCotanWeightsQ},
&IntrinsicGeometryInterface::computeEdgeCotanWeights); addQuantity(vertexAngleDefectsQ, {&halfedgeOppositeAnglesQ},
&IntrinsicGeometryInterface::computeVertexAngleDefects);

  // == Vector fields, angles, and transport
  addQuantity(halfedgeFaceCoordsQ,              {&halfedgeOppositeAnglesQ, &edgeLengthsQ},
&IntrinsicGeometryInterface::computeHalfedgeFaceCoords); addQuantity(faceTransportCoefsQ, {&halfedgeFaceCoordsQ},
&IntrinsicGeometryInterface::computeFaceTransportCoefs); addQuantity(halfedgeOppositeAnglesQ,          {&edgeLengthsQ},
&IntrinsicGeometryInterface::computeHalfedgeOppositeAngles); addQuantity(halfedgeRescaledOppositeAnglesQ,
{&vertexAngleDefectsQ, &halfedgeOppositeAnglesQ},   &IntrinsicGeometryInterface::computeHalfedgeRescaledOppositeAngles);
  addQuantity(halfedgeVertexCoordsQ,            {&halfedgeRescaledOppositeAnglesQ},
&IntrinsicGeometryInterface::computeHalfedgeVertexCoords); addQuantity(vertexTransportCoefsQ, {&halfedgeVertexCoordsQ},
&IntrinsicGeometryInterface::computeVertexTransportCoefs);

  // == Indices
  addQuantity(vertexIndicesQ,             {}, &IntrinsicGeometryInterface::computeVertexIndices);
  addQuantity(interiorVertexIndicesQ,     {}, &IntrinsicGeometryInterface::computeInteriorVertexIndices);
  addQuantity(faceIndicesQ,               {}, &IntrinsicGeometryInterface::computeFaceIndices);
  addQuantity(edgeIndicesQ,               {}, &IntrinsicGeometryInterface::computeEdgeIndices);
  addQuantity(halfedgeIndicesQ,           {}, &IntrinsicGeometryInterface::computeHalfedgeIndices);

  // == Operators
  addQuantity(basicDECOperatorsQ,         {&vertexDualAreasQ, &edgeCotanWeightsQ, &faceAreasQ, &vertexIndicesQ,
&faceIndicesQ, &edgeIndicesQ},      &IntrinsicGeometryInterface::computeBasicDECOperators);
  addQuantity(zeroFormWeakLaplacianQ,     {&basicDECOperatorsQ},
&IntrinsicGeometryInterface::computeZeroFormWeakLaplacian);
  // clang-format on
}
*/


// === Quantity implementations

// Edge lengths
void IntrinsicGeometryInterface::requireEdgeLengths() { edgeLengthsQ.require(); }
void IntrinsicGeometryInterface::unrequireEdgeLengths() { edgeLengthsQ.unrequire(); }

// Face areas
void IntrinsicGeometryInterface::computeFaceAreas() {
  edgeLengthsQ.ensureHave();

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
void IntrinsicGeometryInterface::requireFaceAreas() { faceAreasQ.require(); }
void IntrinsicGeometryInterface::unrequireFaceAreas() { faceAreasQ.unrequire(); }


// Vertex dual area
void IntrinsicGeometryInterface::computeVertexDualAreas() {
  faceAreasQ.ensureHave();

  vertexDualAreas = VertexData<double>(mesh, 0.);

  for (Face f : mesh.faces()) {
    double A = faceAreas[f];
    for (Vertex v : f.adjacentVertices()) {
      vertexDualAreas[v] = A / 3.0;
    }
  }
}
void IntrinsicGeometryInterface::requireVertexDualAreas() { vertexDualAreasQ.require(); }
void IntrinsicGeometryInterface::unrequireVertexDualAreas() { vertexDualAreasQ.unrequire(); }


// Corner angles
void IntrinsicGeometryInterface::computeCornerAngles() {
  edgeLengthsQ.ensureHave();

  cornerAngles = CornerData<double>(mesh);

  for (Corner c : mesh.corners()) {
    Halfedge heA = c.halfedge();
    Halfedge heOpp = heA.next();
    Halfedge heB = heOpp.next();

    GC_SAFETY_ASSERT(heB.next() == heA, "faces mush be triangular");

    double lOpp = edgeLengths[heOpp.edge()];
    double lA = edgeLengths[heA.edge()];
    double lB = edgeLengths[heB.edge()];

    double q = (lA * lA + lB * lB - lOpp * lOpp) / (2. * lA * lB);
    q = clamp(q, -1.0, 1.0);
    double angle = std::acos(q);

    cornerAngles[c] = angle;
  }
}
void IntrinsicGeometryInterface::requireCornerAngles() { cornerAnglesQ.require(); }
void IntrinsicGeometryInterface::unrequireCornerAngles() { cornerAnglesQ.unrequire(); }


// Vertex angle sums
void IntrinsicGeometryInterface::computeVertexAngleSums() {
  cornerAnglesQ.ensureHave();

  vertexAngleSums = VertexData<double>(mesh);
  for (Corner c : mesh.corners()) {
    vertexAngleSums[c.vertex()] += cornerAngles[c];
  }
}
void IntrinsicGeometryInterface::requireVertexAngleSums() { vertexAngleSumsQ.require(); }
void IntrinsicGeometryInterface::unrequireVertexAngleSums() { vertexAngleSumsQ.unrequire(); }


// Corner scaled angles
void IntrinsicGeometryInterface::computeCornerScaledAngles() {
  cornerAnglesQ.ensureHave();
  vertexAngleSumsQ.ensureHave();

  cornerScaledAngles = CornerData<double>(mesh);

  for (Corner c : mesh.corners()) {
    double s = 2.0 * PI / vertexAngleSums[c.vertex()];
    cornerScaledAngles[c] = s * cornerAngles[c];
  }
}
void IntrinsicGeometryInterface::requireCornerScaledAngles() { cornerScaledAnglesQ.require(); }
void IntrinsicGeometryInterface::unrequireCornerScaledAngles() { cornerScaledAnglesQ.unrequire(); }


// Vertex gaussian curvatures
void IntrinsicGeometryInterface::computeVertexGaussianCurvatures() {
  vertexAngleSumsQ.ensureHave();

  vertexGaussianCurvatures = VertexData<double>(mesh, 0);

  for (Vertex v : mesh.vertices()) {
    if (!v.isBoundary()) {
      vertexGaussianCurvatures[v] = 2. * PI - vertexAngleSums[v];
    }
  }
}
void IntrinsicGeometryInterface::requireVertexGaussianCurvatures() { vertexGaussianCurvaturesQ.require(); }
void IntrinsicGeometryInterface::unrequireVertexGaussianCurvatures() { vertexGaussianCurvaturesQ.unrequire(); }

// Face gaussian curvatures
void IntrinsicGeometryInterface::computeFaceGaussianCurvatures() {
  cornerScaledAnglesQ.ensureHave();

  faceGaussianCurvatures = FaceData<double>(mesh);

  for (Face f : mesh.faces()) {

    double angleDefect = -PI;
    Halfedge he = f.halfedge();
    for (int i = 0; i < 3; i++) {
      angleDefect += cornerScaledAngles[he.corner()];
      he = he.next();
    }
    GC_SAFETY_ASSERT(he == f.halfedge(), "faces mush be triangular");

    faceGaussianCurvatures[f] = angleDefect;
  }
}
void IntrinsicGeometryInterface::requireFaceGaussianCurvatures() { faceGaussianCurvaturesQ.require(); }
void IntrinsicGeometryInterface::unrequireFaceGaussianCurvatures() { faceGaussianCurvaturesQ.unrequire(); }

// Halfedge cotan weights
void IntrinsicGeometryInterface::computeHalfedgeCotanWeights() {
  edgeLengthsQ.ensureHave();
  faceAreasQ.ensureHave();

  halfedgeCotanWeights = HalfedgeData<double>(mesh, 0.);

  for (Halfedge he : mesh.interiorHalfedges()) {

    Halfedge heF = he;
    double l_ij = edgeLengths[heF.edge()];
    heF = heF.next();
    double l_jk = edgeLengths[heF.edge()];
    heF = heF.next();
    double l_ki = edgeLengths[heF.edge()];
    heF = heF.next();

    GC_SAFETY_ASSERT(heF == he, "faces mush be triangular");

    double area = faceAreas[he.face()];
    double cotValue = (-l_ij * l_ij + l_jk * l_jk + l_ki * l_ki) / (4. * area);
    halfedgeCotanWeights[he] = cotValue / 2;
  }
}
void IntrinsicGeometryInterface::requireHalfedgeCotanWeights() { halfedgeCotanWeightsQ.require(); }
void IntrinsicGeometryInterface::unrequireHalfedgeCotanWeights() { halfedgeCotanWeightsQ.unrequire(); }

// Edge cotan weights
void IntrinsicGeometryInterface::computeEdgeCotanWeights() {
  edgeLengthsQ.ensureHave();
  faceAreasQ.ensureHave();

  edgeCotanWeights = EdgeData<double>(mesh, 0.);

  for (Edge e : mesh.edges()) {
    double cotSum = 0.;

    { // First halfedge-- always real
      Halfedge he = e.halfedge();
      double l_ij = edgeLengths[he.edge()];
      he = he.next();
      double l_jk = edgeLengths[he.edge()];
      he = he.next();
      double l_ki = edgeLengths[he.edge()];
      he = he.next();
      GC_SAFETY_ASSERT(he == e.halfedge(), "faces mush be triangular");
      double area = faceAreas[he.face()];
      double cotValue = (-l_ij * l_ij + l_jk * l_jk + l_ki * l_ki) / (4. * area);
      cotSum += cotValue / 2;
    }

    if (e.halfedge().twin().isInterior()) { // Second halfedge
      Halfedge he = e.halfedge().twin();
      double l_ij = edgeLengths[he.edge()];
      he = he.next();
      double l_jk = edgeLengths[he.edge()];
      he = he.next();
      double l_ki = edgeLengths[he.edge()];
      double area = faceAreas[he.face()];
      double cotValue = (-l_ij * l_ij + l_jk * l_jk + l_ki * l_ki) / (4. * area);
      cotSum += cotValue / 2;
    }

    edgeCotanWeights[e] = cotSum;
  }
}
void IntrinsicGeometryInterface::requireEdgeCotanWeights() { edgeCotanWeightsQ.require(); }
void IntrinsicGeometryInterface::unrequireEdgeCotanWeights() { edgeCotanWeightsQ.unrequire(); }

/*
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


void IntrinsicGeometryInterface::computeZeroFormWeakLaplacian() {
  zeroFormWeakLaplacian = d0.transpose() * hodge1 * d0;
}
*/


} // namespace surface
} // namespace geometrycentral
