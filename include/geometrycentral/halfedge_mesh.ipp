#pragma once

namespace geometrycentral {

// Methods for getting number of mesh elements

inline size_t HalfedgeMesh::nHalfedges(void) const { return nRealHalfedgesCount; }
inline size_t HalfedgeMesh::nCorners(void) const { return nRealHalfedgesCount; }
inline size_t HalfedgeMesh::nVertices(void) const { return nVerticesCount; }
inline size_t HalfedgeMesh::nEdges(void) const { return nEdgesCount; }
inline size_t HalfedgeMesh::nFaces(void) const { return nFacesCount; }
inline size_t HalfedgeMesh::nBoundaryLoops(void) const { return rawBoundaryLoops.size(); }
inline size_t HalfedgeMesh::nImaginaryHalfedges(void) const { return nImaginaryHalfedgesCount; }

inline size_t HalfedgeMesh::nHalfedgesCapacity(void) const { return rawHalfedges.capacity(); }
inline size_t HalfedgeMesh::nVerticesCapacity(void) const { return rawVertices.capacity(); }
inline size_t HalfedgeMesh::nEdgesCapacity(void) const { return rawEdges.capacity(); }
inline size_t HalfedgeMesh::nFacesCapacity(void) const { return rawFaces.capacity(); }

// Methods for iterating over mesh elements w/ range-based for loops ===========

inline HalfedgePtrSet HalfedgeMesh::realHalfedges(void) {
  size_t nH = rawHalfedges.size();
  HalfedgePtr beginptr(&rawHalfedges[0]);
  HalfedgePtr endptr(&rawHalfedges[nH]);

  return HalfedgePtrSet(beginptr, endptr, HalfedgeSetType::Real);
}

inline HalfedgePtrSet HalfedgeMesh::imaginaryHalfedges(void) {
  size_t nH = rawHalfedges.size();
  HalfedgePtr beginptr(&rawHalfedges[0]);
  HalfedgePtr endptr(&rawHalfedges[nH]);

  return HalfedgePtrSet(beginptr, endptr, HalfedgeSetType::Imaginary);
}

inline HalfedgePtrSet HalfedgeMesh::allHalfedges(void) {
  size_t nH = rawHalfedges.size();
  HalfedgePtr beginptr(&rawHalfedges[0]);
  HalfedgePtr endptr(&rawHalfedges[nH]);

  return HalfedgePtrSet(beginptr, endptr, HalfedgeSetType::All);
}

inline CornerPtrSet HalfedgeMesh::corners(void) {
  size_t nC = rawHalfedges.size();
  CornerPtr beginptr(&rawHalfedges[0]);
  CornerPtr endptr(&rawHalfedges[nC]);

  return CornerPtrSet(beginptr, endptr);
}

inline VertexPtrSet HalfedgeMesh::vertices(void) {
  size_t nV = rawVertices.size();
  VertexPtr beginptr{&rawVertices[0]};
  VertexPtr endptr{&rawVertices[nV]};

  return VertexPtrSet(beginptr, endptr);
}

inline EdgePtrSet HalfedgeMesh::edges(void) {
  size_t nE = rawEdges.size();
  EdgePtr beginptr{&rawEdges[0]};
  EdgePtr endptr{&rawEdges[nE]};

  return EdgePtrSet(beginptr, endptr);
}

inline FacePtrSet HalfedgeMesh::faces(void) {
  size_t nF = rawFaces.size();
  FacePtr beginptr{&rawFaces[0]};
  FacePtr endptr{&rawFaces[nF]};

  return FacePtrSet(beginptr, endptr);
}

inline BoundaryPtrSet HalfedgeMesh::boundaryLoops(void) {
  size_t nBL = rawBoundaryLoops.size();
  BoundaryLoopPtr beginptr{&rawBoundaryLoops[0]};
  BoundaryLoopPtr endptr{&rawBoundaryLoops[nBL]};

  return BoundaryPtrSet(beginptr, endptr);
}

// Methods for accessing elements by index =====================================
// Note that these are only valid when the mesh is compressed.

inline HalfedgePtr HalfedgeMesh::halfedge(size_t index) { return HalfedgePtr{&rawHalfedges[index]}; }

inline CornerPtr HalfedgeMesh::corner(size_t index) { return CornerPtr{&rawHalfedges[index]}; }

inline VertexPtr HalfedgeMesh::vertex(size_t index) { return VertexPtr{&rawVertices[index]}; }

inline EdgePtr HalfedgeMesh::edge(size_t index) { return EdgePtr{&rawEdges[index]}; }

inline FacePtr HalfedgeMesh::face(size_t index) { return FacePtr{&rawFaces[index]}; }

inline BoundaryLoopPtr HalfedgeMesh::boundaryLoop(size_t index) { return BoundaryLoopPtr{&rawBoundaryLoops[index]}; }

// Misc utility methods =====================================

inline bool HalfedgeMesh::isCompressed() { return isCompressedFlag; }
inline bool HalfedgeMesh::isCanonical() { return isCanonicalFlag; }

} // namespace geometrycentral
