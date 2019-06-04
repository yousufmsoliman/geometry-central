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
  inline void requireEdgeLengths();
  EdgeData<double> edgeLengths;

  // Face areas
  inline void requireFaceAreas();
  FaceData<double> faceAreas;

  /*
  // Face areas
  inline void requireFaceAreas() { faceAreasQ.require(); }
  FaceData<double> faceAreas;

  // Vertex dual areas
  inline void requireVertexDualAreas() { vertexDualAreasQ.require(); }
  VertexData<double> vertexDualAreas;

  // Edge lengths
  inline void requireEdgeLengths() { edgeLengthsQ.require(); }
  EdgeData<double> edgeLengths;

  // Halfedge cotan weights
  inline void requireHalfedgeCotanWeights() { halfedgeCotanWeightsQ.require(); }
  HalfedgeData<double> halfedgeCotanWeights;

  // Edge cotan weights
  inline void requireEdgeCotanWeights() { edgeCotanWeightsQ.require(); }
  EdgeData<double> edgeCotanWeights;

  // Angle defect at vertices
  inline void requireVertexAngleDefects() { vertexAngleDefectsQ.require(); }
  VertexData<double> vertexAngleDefects;


  // == Vector fields, angles, and transport

  // The coordinate of each halfedge in the basis of he.face()
  // NOTE: These HAVE magnitude, unlike the vertex version (confusingly)
  inline void requireHalfedgeFaceCoords() { halfedgeFaceCoordsQ.require(); }
  HalfedgeData<Complex> halfedgeFaceCoords;

  // Transport an intrinsic vector field in he.face() to he.twin().face() by multiplying
  // by complex z = e^(theta I) given here
  inline void requireFaceTransportCoefs() { faceTransportCoefsQ.require(); }
  HalfedgeData<Complex> faceTransportCoefs;

  // Halfedge opposite angles
  inline void requireHalfedgeOppositeAngles() { halfedgeOppositeAnglesQ.require(); }
  HalfedgeData<double> halfedgeOppositeAngles;

  // Halfedge opposite angles (scaled by the angle defect to sum to 2 PI at each vertex)
  inline void requireHalfedgeRescaledOppositeAngles() { halfedgeRescaledOppositeAnglesQ.require(); }
  HalfedgeData<double> halfedgeRescaledOppositeAngles;

  // The coordinate of each halfedge in the basis of he.vertex(), rescaled so the sum around each vertex is 2*PI
  inline void requireHalfedgeVertexCoords() { halfedgeVertexCoordsQ.require(); }
  HalfedgeData<Complex> halfedgeVertexCoords;

  // Transport an intrinsic vector field in he.vertex() to he.twin().vertex() by multiplying
  // by complex z = e^(theta I) given here
  inline void requireVertexTransportCoefs() { vertexTransportCoefsQ.require(); }
  HalfedgeData<Complex> vertexTransportCoefs;


  // == Operators
  // Note: These don't quite follow the usual naming scheme, for the sake of grouping common operators
  // TODO factorizations?

  // All of the basic DEC operators and their inverses
  inline void requireBasicDECOperators() { basicDECOperatorsQ.require(); }
  Eigen::SparseMatrix<double> d0, d1, hodge0, hodge1, hodge2;
  Eigen::SparseMatrix<double> hodge0Inv, hodge1Inv, hodge2Inv;

  // Cotan-laplace operator
  // Remember, this DOES NOT include the mass matrix (hodge0)
  inline void requireZeroFormWeakLaplacian() { zeroFormWeakLaplacianQ.require(); }
  Eigen::SparseMatrix<double> zeroFormWeakLaplacian;
  */

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
