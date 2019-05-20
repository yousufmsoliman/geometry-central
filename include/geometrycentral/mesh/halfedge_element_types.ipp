#pragma once

// Implementations for halfedge_mesh_types.ipp

namespace geometrycentral {
namespace halfedge_mesh {

// clang-format off

// ==========================================================
// ================      Base Element      ==================
// ==========================================================

template<typename T> 
inline bool Element<T>::operator==(const Element<T>& other) const { return ind == other.ind; }
template<typename T> 
inline bool Element<T>::operator!=(const Element<T>& other) const { return !(*this == other); }
template<typename T> 
inline bool Element<T>::operator>(const Element<T>& other) const { return ind > other.ind; }
template<typename T> 
inline bool Element<T>::operator>=(const Element<T>& other) const { return ind >= other.ind; }
template<typename T> 
inline bool Element<T>::operator<(const Element<T>& other) const { return ind < other.ind; }
template<typename T> 
inline bool Element<T>::operator<=(const Element<T>& other) const { return ind <= other.ind; }

template <typename T>
size_t Element<T>::getIndex() const { return ind; }

template <typename T>
HalfedgeMesh* Element<T>::getMesh() const { return mesh; }

template <typename T>
inline ::std::ostream& operator<<(::std::ostream& output, const Element<T>& e) {
  output << typeShortName<T>() << "_" << e.ind;
  return output;
}


// Base iterators
template <typename F>
inline RangeIteratorBase<F>::RangeIteratorBase(HalfedgeMesh* mesh_, size_t iStart_, size_t iEnd_) : mesh(mesh_), iCurr(iStart_), iEnd(iEnd_) {
  if (iCurr != iEnd && F::elementOkay(*mesh, iCurr)) {
    this->operator++();
  }
}

template <typename F>
inline const RangeIteratorBase<F>& RangeIteratorBase<F>::operator++() {
  iCurr++;
  while (iCurr != iEnd && F::elementOkay(*mesh, iCurr)) {
    iCurr++;
  }
  return *this;
}

template <typename F>
inline bool RangeIteratorBase<F>::operator==(const RangeIteratorBase<F>& other) const {
	return iCurr == other.iCurr;
}

template <typename F>
inline bool RangeIteratorBase<F>::operator!=(const RangeIteratorBase<F>& other) const {
	return !(*this == other);
}

template <typename F>
inline typename F::Etype RangeIteratorBase<F>::operator*() const { return typename F::Etype(mesh, iCurr); }

template <typename F>
RangeSetBase<F>::RangeSetBase(HalfedgeMesh* mesh_, size_t iStart_, size_t iEnd_) : mesh(mesh_), iStart(iStart_), iEnd(iEnd_) {}  

template <typename F>
inline RangeIteratorBase<F> RangeSetBase<F>::begin() const { return RangeIteratorBase<F>(mesh, iStart, iEnd); }

template <typename F>
inline RangeIteratorBase<F> RangeSetBase<F>::end() const { return RangeIteratorBase<F>(mesh, iEnd, iEnd); }

// ==========================================================
// ================        Vertex          ==================
// ==========================================================

// Navigators
inline Halfedge Vertex::halfedge() const    { return Halfedge(mesh, mesh->vHalfedge[ind]); }
inline Corner Vertex::corner() const        { return halfedge().corner(); }

// Properties
inline bool Vertex::isBoundary() const { return halfedge().twin().isInterior(); }

// Iterators
inline VertexIncomingHalfedgeSet Vertex::incomingHalfedges() const { return VertexIncomingHalfedgeSet(halfedge().twin()); }
inline VertexOutgoingHalfedgeSet Vertex::outgoingHalfedges() const { return VertexOutgoingHalfedgeSet(halfedge()); }
inline VertexAdjacentVertexSet Vertex::adjacentVertices() const { return VertexAdjacentVertexSet(halfedge().twin()); }
inline VertexAdjacentFaceSet Vertex::adjacentFaces() const { return VertexAdjacentFaceSet(halfedge()); }
inline VertexAdjacentEdgeSet Vertex::adjacentEdges() const { return VertexAdjacentEdgeSet(halfedge()); }
inline VertexAdjacentCornerSet Vertex::adjacentCorners() const {
  Halfedge h = halfedge();
  if (!h.isInterior()) h = h.twin().next();
  return VertexAdjacentCornerSet(h);
}

// Range iterators
inline bool VertexRangeF::elementOkay(const HalfedgeMesh& mesh, size_t ind) {
  return !mesh.vertexIsDead(ind);
}

// ==========================================================
// ================        Halfedge        ==================
// ==========================================================

// Navigators
inline Halfedge Halfedge::twin() const  { return Halfedge(mesh, HalfedgeMesh::heTwin(ind)); }
inline Halfedge Halfedge::next() const  { return Halfedge(mesh, mesh->heNext[ind]); }
inline Vertex Halfedge::vertex() const  { return Vertex(mesh, mesh->heVertex[ind]); }
inline Edge Halfedge::edge() const      { return Edge(mesh, HalfedgeMesh::heEdge(ind)); }
inline Face Halfedge::face() const      { return Face(mesh, mesh->heFace[ind]); }
inline Corner Halfedge::corner() const  { return Corner(mesh, ind); }

// Properties
inline bool Halfedge::isInterior() const { return  mesh->heIsInterior(ind); }

// Range iterators
inline bool HalfedgeRangeF::elementOkay(const HalfedgeMesh& mesh, size_t ind) {
  return !mesh.halfedgeIsDead(ind);
}
inline bool HalfedgeInteriorRangeF::elementOkay(const HalfedgeMesh& mesh, size_t ind) {
  return !mesh.halfedgeIsDead(ind) && mesh.heIsInterior(ind);
}
inline bool HalfedgeExteriorRangeF::elementOkay(const HalfedgeMesh& mesh, size_t ind) {
  return !mesh.halfedgeIsDead(ind) && !mesh.heIsInterior(ind);
}

// ==========================================================
// ================        Corner          ==================
// ==========================================================

// Navigators
inline Halfedge Corner::halfedge() const { return Halfedge(mesh, ind); }
inline Vertex Corner::vertex() const { return halfedge().vertex(); }
inline Face Corner::face() const { return halfedge().face(); }

// Range iterators
inline bool CornerRangeF::elementOkay(const HalfedgeMesh& mesh, size_t ind) {
  return !mesh.halfedgeIsDead(ind) && mesh.heIsInterior(ind);
}

// ==========================================================
// ================          Edge          ==================
// ==========================================================

// Navigators
inline Halfedge Edge::halfedge() const { return Halfedge(mesh, HalfedgeMesh::eHalfedge(ind)); }

// Properties
inline bool Edge::isBoundary() const { return !halfedge().isInterior() || !halfedge().twin().isInterior(); }

// Range iterators
inline bool EdgeRangeF::elementOkay(const HalfedgeMesh& mesh, size_t ind) {
  return !mesh.edgeIsDead(ind) && mesh.heIsInterior(ind);
}


// ==========================================================
// ================          Face          ==================
// ==========================================================

// Navigation
inline Halfedge Face::halfedge() const { return Halfedge(mesh, mesh->fHalfedge[ind]); }
inline BoundaryLoop Face::asBoundaryLoop() const { 
  GC_SAFETY_ASSERT(isBoundaryLoop(), "face must be boundary loop to call asBoundaryLoop()")
  return BoundaryLoop(mesh, mesh->faceIndToBoundaryLoopInd(ind)); 
}

// Properties
inline size_t Face::degree() const {
  size_t k = 0;
  for (Halfedge h : adjacentHalfedges()) { k++; }
  return k;
}

inline bool Face::isBoundaryLoop() const { return mesh->faceIsBoundaryLoop(ind); }

// Iterators
inline FaceAdjacentHalfedgeSet Face::adjacentHalfedges() const { return FaceAdjacentHalfedgeSet(halfedge()); }
inline FaceAdjacentVertexSet Face::adjacentVertices() const { return FaceAdjacentVertexSet(halfedge()); }
inline FaceAdjacentFaceSet Face::adjacentFaces() const { return FaceAdjacentFaceSet(halfedge()); }
inline FaceAdjacentEdgeSet Face::adjacentEdges() const { return FaceAdjacentEdgeSet(halfedge()); }
inline FaceAdjacentCornerSet Face::adjacentCorners() const { return FaceAdjacentCornerSet(halfedge()); }


// Range iterators
inline bool FaceRangeF::elementOkay(const HalfedgeMesh& mesh, size_t ind) {
  return !mesh.faceIsDead(ind);
}

// ==========================================================
// ================     Boundary Loop      ==================
// ==========================================================

inline Halfedge BoundaryLoop::halfedge() const { return asFace().halfedge(); }
inline Face BoundaryLoop::asFace() const { return Face(mesh, mesh->boundaryLoopIndToFaceInd(ind)); }
inline size_t BoundaryLoop::degree() const {
  size_t k = 0;
  for (Halfedge h : adjacentHalfedges()) { k++; }
  return k;
}

inline BoundaryLoopAdjacentHalfedgeSet BoundaryLoop::adjacentHalfedges() const { return BoundaryLoopAdjacentHalfedgeSet(halfedge()); }
inline BoundaryLoopAdjacentVertexSet BoundaryLoop::adjacentVertices() const { return BoundaryLoopAdjacentVertexSet(halfedge()); }
inline BoundaryLoopAdjacentEdgeSet BoundaryLoop::adjacentEdges() const { return BoundaryLoopAdjacentEdgeSet(halfedge()); }

// Range iterators
inline bool BoundaryLoopRangeF::elementOkay(const HalfedgeMesh& mesh, size_t ind) {
  return !mesh.faceIsDead(mesh.boundaryLoopIndToFaceInd(ind));
}


// clang-format on

} // namespace halfedge_mesh
} // namespace geometrycentral

namespace std {

template <typename T>
struct hash<geometrycentral::halfedge_mesh::Element<T>> {
  std::size_t operator()(const geometrycentral::halfedge_mesh::Element<T>& e) const {
    return std::hash<size_t>{}(e.getIndex());
  }
};
// template <>
// struct hash<geometrycentral::Vertex> {
// std::size_t operator()(const geometrycentral::Vertex& v) const { return std::hash<size_t>{}(v.ptr->ID); }
//};

} // namespace std
