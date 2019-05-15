#pragma once

namespace geometrycentral {

// Methods for getting number of mesh elements

inline size_t HalfedgeMesh::nRealHalfedges(void) const { return nRealHalfedgesCount; }
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

inline HalfedgeSet HalfedgeMesh::realHalfedges(void) {
  size_t nH = rawHalfedges.size();
  Halfedge beginptr(&rawHalfedges[0]);
  Halfedge endptr(&rawHalfedges[nH]);

  return HalfedgeSet(beginptr, endptr, HalfedgeSetType::Real);
}

inline HalfedgeSet HalfedgeMesh::imaginaryHalfedges(void) {
  size_t nH = rawHalfedges.size();
  Halfedge beginptr(&rawHalfedges[0]);
  Halfedge endptr(&rawHalfedges[nH]);

  return HalfedgeSet(beginptr, endptr, HalfedgeSetType::Imaginary);
}

inline HalfedgeSet HalfedgeMesh::allHalfedges(void) {
  size_t nH = rawHalfedges.size();
  Halfedge beginptr(&rawHalfedges[0]);
  Halfedge endptr(&rawHalfedges[nH]);

  return HalfedgeSet(beginptr, endptr, HalfedgeSetType::All);
}

inline CornerSet HalfedgeMesh::corners(void) {
  size_t nC = rawHalfedges.size();
  Corner beginptr(&rawHalfedges[0]);
  Corner endptr(&rawHalfedges[nC]);

  return CornerSet(beginptr, endptr);
}

inline VertexSet HalfedgeMesh::vertices(void) {
  size_t nV = rawVertices.size();
  Vertex beginptr{&rawVertices[0]};
  Vertex endptr{&rawVertices[nV]};

  return VertexSet(beginptr, endptr);
}

inline EdgeSet HalfedgeMesh::edges(void) {
  size_t nE = rawEdges.size();
  Edge beginptr{&rawEdges[0]};
  Edge endptr{&rawEdges[nE]};

  return EdgeSet(beginptr, endptr);
}

inline FaceSet HalfedgeMesh::faces(void) {
  size_t nF = rawFaces.size();
  Face beginptr{&rawFaces[0]};
  Face endptr{&rawFaces[nF]};

  return FaceSet(beginptr, endptr);
}

inline BoundarySet HalfedgeMesh::boundaryLoops(void) {
  size_t nBL = rawBoundaryLoops.size();
  BoundaryLoop beginptr{&rawBoundaryLoops[0]};
  BoundaryLoop endptr{&rawBoundaryLoops[nBL]};

  return BoundarySet(beginptr, endptr);
}

// Methods for accessing elements by index =====================================
// Note that these are only valid when the mesh is compressed.

inline Halfedge HalfedgeMesh::halfedge(size_t index) { return Halfedge{&rawHalfedges[index]}; }

inline Corner HalfedgeMesh::corner(size_t index) { return Corner{&rawHalfedges[index]}; }

inline Vertex HalfedgeMesh::vertex(size_t index) { return Vertex{&rawVertices[index]}; }

inline Edge HalfedgeMesh::edge(size_t index) { return Edge{&rawEdges[index]}; }

inline Face HalfedgeMesh::face(size_t index) { return Face{&rawFaces[index]}; }

inline BoundaryLoop HalfedgeMesh::boundaryLoop(size_t index) { return BoundaryLoop{&rawBoundaryLoops[index]}; }

// Misc utility methods =====================================

inline bool HalfedgeMesh::isCompressed() { return isCompressedFlag; }
inline bool HalfedgeMesh::isCanonical() { return isCanonicalFlag; }

} // namespace geometrycentral
