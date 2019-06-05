#pragma once

#include "geometrycentral/surface/base_geometry_interface.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/utilities/vector2.h"

#include <Eigen/SparseCore>

namespace geometrycentral {
namespace surface {


class IntrinsicGeometryInterface : public BaseGeometryInterface {

protected:
  // Constructor is protected, because this class is an interface which is not meant to be instantiated directly.
  // Instantiate it via some realization which encapsulates input data, like EdgeLengthGeometry or
  // VertexPositionGeometry.
  IntrinsicGeometryInterface(HalfedgeMesh& mesh_);
  virtual ~IntrinsicGeometryInterface() {}

public:
  // == Lengths, areas, and angles

  // Edge lengths
  EdgeData<double> edgeLengths;
  void requireEdgeLengths();
  void unrequireEdgeLengths();

  // Face areas
  FaceData<double> faceAreas;
  void requireFaceAreas();
  void unrequireFaceAreas();

  // Vertex dual areas
  VertexData<double> vertexDualAreas;
  void requireVertexDualAreas();
  void unrequireVertexDualAreas();

  // Corner angles
  CornerData<double> cornerAngles;
  void requireCornerAngles();
  void unrequireCornerAngles();

  // Vertex angle sums
  VertexData<double> vertexAngleSums;
  void requireVertexAngleSums();
  void unrequireVertexAngleSums();

  // Corner scaled angles
  CornerData<double> cornerScaledAngles;
  void requireCornerScaledAngles();
  void unrequireCornerScaledAngles();

  // Vertex gaussian curvature
  VertexData<double> vertexGaussianCurvatures;
  void requireVertexGaussianCurvatures();
  void unrequireVertexGaussianCurvatures();

  // Face gaussian curvature
  FaceData<double> faceGaussianCurvatures;
  void requireFaceGaussianCurvatures();
  void unrequireFaceGaussianCurvatures();

  // Halfedge cotan weight
  HalfedgeData<double> halfedgeCotanWeights;
  void requireHalfedgeCotanWeights();
  void unrequireHalfedgeCotanWeights();

  // Edge cotan weight
  EdgeData<double> edgeCotanWeights;
  void requireEdgeCotanWeights();
  void unrequireEdgeCotanWeights();


  // == Tangent vectors and transport
 
  // Halfedge vectors in face tangent space
  HalfedgeData<Vector2> halfedgeVectorsInFace;
  void requireHalfedgeVectorsInFace();
  void unrequireHalfedgeVectorsInFace();

  // Face tangent vector transport across halfedges
  HalfedgeData<Vector2> transportVectorsAcrossHalfedge;
  void requireTransportVectorsAcrossHalfedge();
  void unrequireTransportVectorsAcrossHalfedge();
  
  // Halfedge vectors in vertex tangent space
  HalfedgeData<Vector2> halfedgeVectorsInVertex;
  void requireHalfedgeVectorsInVertex();
  void unrequireHalfedgeVectorsInVertex();

  // Vertex transport across halfedges
  HalfedgeData<Vector2> transportVectorsAlongHalfedge;
  void requireTransportVectorsAlongHalfedge();
  void unrequireTransportVectorsAlongHalfedge();

protected:

  // == Lengths, areas, and angles

  // Edge lengths
  // Note that computeEdgeLengths() is pure virtual: some input data class which extends this interface must supply a
  // method for computing edge lengths (EdgeLengthGeometry serves this purpose)
  DependentQuantityD<EdgeData<double>> edgeLengthsQ;
  virtual void computeEdgeLengths() = 0;

  // Face areas
  DependentQuantityD<FaceData<double>> faceAreasQ;
  virtual void computeFaceAreas();

  // Vertex dual area
  DependentQuantityD<VertexData<double>> vertexDualAreasQ;
  virtual void computeVertexDualAreas();
  
  // Corner angles
  DependentQuantityD<CornerData<double>> cornerAnglesQ;
  virtual void computeCornerAngles();

  // Vertex angle sums
  DependentQuantityD<VertexData<double>> vertexAngleSumsQ;
  virtual void computeVertexAngleSums();  
  
  // Corner scaled angles
  DependentQuantityD<CornerData<double>> cornerScaledAnglesQ;
  virtual void computeCornerScaledAngles();
   
  // Vertex gaussian curvature
  DependentQuantityD<VertexData<double>> vertexGaussianCurvaturesQ;
  virtual void computeVertexGaussianCurvatures();
  
  // Face gaussian curvature
  DependentQuantityD<FaceData<double>> faceGaussianCurvaturesQ;
  virtual void computeFaceGaussianCurvatures();
  
  // Halfedge cotan weight
  DependentQuantityD<HalfedgeData<double>> halfedgeCotanWeightsQ;
  virtual void computeHalfedgeCotanWeights();
  
  // Edge cotan weight
  DependentQuantityD<EdgeData<double>> edgeCotanWeightsQ;
  virtual void computeEdgeCotanWeights();


  // == Tangent vectors and transport
  
  // Halfedge vectors in face
  DependentQuantityD<HalfedgeData<Vector2>> halfedgeVectorsInFaceQ;
  virtual void computeHalfedgeVectorsInFace();
  
  // Face tangent vector transport across halfedges
  DependentQuantityD<HalfedgeData<Vector2>> transportVectorsAcrossHalfedgeQ;
  virtual void computeTransportVectorsAcrossHalfedge();
  
  // Halfedge vectors in vertex tangent space
  DependentQuantityD<HalfedgeData<Vector2>> halfedgeVectorsInVertexQ;
  virtual void computeHalfedgeVectorsInVertex();
  
  // Vertex transport across halfedges
  DependentQuantityD<HalfedgeData<Vector2>> transportVectorsAlongHalfedgeQ;
  virtual void computeTransportVectorsAlongHalfedge();

  /*
  // == Basic geometric quantities

  DependentQuantity faceAreasQ;
  virtual void computeFaceAreas();

  DependentQuantity vertexDualAreasQ;
  virtual void computeVertexDualAreas();

  DependentQuantity edgeLengthsQ;
  virtual void computeEdgeLengths() = 0;

  DependentQuantity halfedgeCotanWeightsQ;
  virtual void computeHalfedgeCotanWeights();

  DependentQuantity edgeCotanWeightsQ;
  virtual void computeEdgeCotanWeights();

  DependentQuantity vertexAngleDefectsQ;
  virtual void computeVertexAngleDefects();

  // == Vector fields, angles, and transport

  DependentQuantity halfedgeFaceCoordsQ;
  virtual void computeHalfedgeFaceCoords();

  DependentQuantity faceTransportCoefsQ;
  virtual void computeFaceTransportCoefs();

  DependentQuantity halfedgeOppositeAnglesQ;
  virtual void computeHalfedgeOppositeAngles();

  DependentQuantity halfedgeRescaledOppositeAnglesQ;
  virtual void computeHalfedgeRescaledOppositeAngles();

  DependentQuantity halfedgeVertexCoordsQ;
  virtual void computeHalfedgeVertexCoords();

  DependentQuantity vertexTransportCoefsQ;
  virtual void computeVertexTransportCoefs();


  // == Indices

  DependentQuantity vertexIndicesQ;
  virtual void computeVertexIndices();

  DependentQuantity interiorVertexIndicesQ;
  virtual void computeInteriorVertexIndices();

  DependentQuantity faceIndicesQ;
  virtual void computeFaceIndices();

  DependentQuantity edgeIndicesQ;
  virtual void computeEdgeIndices();

  DependentQuantity halfedgeIndicesQ;
  virtual void computeHalfedgeIndices();

  // == Operators

  DependentQuantity basicDECOperatorsQ;
  virtual void computeBasicDECOperators();

  DependentQuantity zeroFormWeakLaplacianQ;
  virtual void computeZeroFormWeakLaplacian();

  */
};

} // namespace surface
} // namespace geometrycentral
