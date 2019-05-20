#pragma once

namespace geometrycentral {
namespace halfedge_mesh {

// ==========================================================
// ================     Base  Iterator     ==================
// ==========================================================

// Note: the class below may seem a bit weird, in that we call advance() lazily in operator* rather than proactively in
// operator++ (and occaisonally in operator++ in case operator* was never called). We do this to avoid repeatedly
// advancing to hte first valid element on std::end(set). The intuition is that the iterators will still be efficient
// even if advancing to the first valid element is expensive (though we do assume that calling isValid() is cheap).
//
// If we one day really want to optimize performance, it might be a good idea to benchmark these to figure out how much
// it costs to satisfy isValid(), which isn't actually used by many iterators in practice.

template <typename N>
inline NavigationIteratorBase<N>::NavigationIteratorBase(HalfedgeMesh* mesh_, typename N::Etype e, bool justStarted_)
    : mesh(mesh_), state(e), justStarted(justStarted_) {}

template <typename N>
inline const NavigationIteratorBase<N>& NavigationIteratorBase<N>::operator++() {
  while (!state.isValid()) { // in case operator* never got called
    state.advance();
  }
  state.advance();
  justStarted = false;
  return *this;
}

template <typename N>
inline bool NavigationIteratorBase<N>::operator==(const NavigationIteratorBase<N>& other) const {
  return justStarted == other.justStarted && state.currE == other.state.currE;
}

template <typename N>
inline bool NavigationIteratorBase<N>::operator!=(const NavigationIteratorBase<N>& other) const {
  return !(*this == other);
}

template <typename N>
inline typename N::Rtype NavigationIteratorBase<N>::operator*() const {
  while (!state.isValid()) {
    state.advance();
  }
  return state.getCurrent();
}

template <typename N>
NavigationSetBase<N>::NavigationSetBase(HalfedgeMesh* mesh_, typename N::Etype firstE_)
    : mesh(mesh_), firstE(firstE_) {}

template <typename N>
inline NavigationIteratorBase<N> NavigationSetBase<N>::begin() const {
  return NavigationIteratorBase<N>(mesh, firstE, true);
}

template <typename N>
inline NavigationIteratorBase<N> NavigationSetBase<N>::end() const {
  return NavigationIteratorBase<N>(mesh, firstE, false);
}


// ==========================================================
// ================    Vertex Iterators    ==================
// ==========================================================

inline void VertexAdjacentVertexNavigator::advance() { currE = currE.twin().next(); }
inline bool VertexAdjacentVertexNavigator::isValid() const { return true; }
inline Vertex VertexAdjacentVertexNavigator::getCurrent() const { return currE.twin().vertex(); }

inline void VertexIncomingHalfedgeNavigator::advance() { currE = currE.twin().next(); }
inline bool VertexIncomingHalfedgeNavigator::isValid() const { return true; }
inline Halfedge VertexIncomingHalfedgeNavigator::getCurrent() const { return currE.twin(); }

inline void VertexOutgoingHalfedgeNavigator::advance() { currE = currE.twin().next(); }
inline bool VertexOutgoingHalfedgeNavigator::isValid() const { return true; }
inline Halfedge VertexOutgoingHalfedgeNavigator::getCurrent() const { return currE; }

inline void VertexAdjacentCornerNavigator::advance() { currE = currE.twin().next(); }
inline bool VertexAdjacentCornerNavigator::isValid() const { return true; }
inline Corner VertexAdjacentCornerNavigator::getCurrent() const { return currE.corner(); }

inline void VertexAdjacentEdgeNavigator::advance() { currE = currE.twin().next(); }
inline bool VertexAdjacentEdgeNavigator::isValid() const { return true; }
inline Edge VertexAdjacentEdgeNavigator::getCurrent() const { return currE.edge(); }

inline void VertexAdjacentFaceNavigator::advance() { currE = currE.twin().next(); }
inline bool VertexAdjacentFaceNavigator::isValid() const { return currE.isInterior(); }
inline Face VertexAdjacentFaceNavigator::getCurrent() const { return currE.face(); }


// ==========================================================
// ================     Face Iterators     ==================
// ==========================================================

inline void FaceAdjacentVertexNavigator::advance() { currE = currE.next(); }
inline bool FaceAdjacentVertexNavigator::isValid() const { return true; }
inline Vertex FaceAdjacentVertexNavigator::getCurrent() const { return currE.vertex(); }

inline void FaceAdjacentHalfedgeNavigator::advance() { currE = currE.next(); }
inline bool FaceAdjacentHalfedgeNavigator::isValid() const { return true; }
inline Halfedge FaceAdjacentHalfedgeNavigator::getCurrent() const { return currE; }

inline void FaceAdjacentCornerNavigator::advance() { currE = currE.next(); }
inline bool FaceAdjacentCornerNavigator::isValid() const { return true; }
inline Corner FaceAdjacentCornerNavigator::getCurrent() const { return currE.corner(); }

inline void FaceAdjacentEdgeNavigator::advance() { currE = currE.next(); }
inline bool FaceAdjacentEdgeNavigator::isValid() const { return true; }
inline Edge FaceAdjacentEdgeNavigator::getCurrent() const { return currE.edge(); }

// TODO I think this has a problem and will loop forever on a single isolated face
inline void FaceAdjacentFaceNavigator::advance() { currE = currE.next(); }
inline bool FaceAdjacentFaceNavigator::isValid() const { return currE.twin().isInterior(); }
inline Face FaceAdjacentFaceNavigator::getCurrent() const { return currE.twin().face(); }

// ==========================================================
// ==============   Boundary Loop Iterators   ===============
// ==========================================================

inline void BoundaryLoopAdjacentVertexNavigator::advance() { currE = currE.next(); }
inline bool BoundaryLoopAdjacentVertexNavigator::isValid() const { return true; }
inline Vertex BoundaryLoopAdjacentVertexNavigator::getCurrent() const { return currE.vertex(); }

inline void BoundaryLoopAdjacentHalfedgeNavigator::advance() { currE = currE.next(); }
inline bool BoundaryLoopAdjacentHalfedgeNavigator::isValid() const { return true; }
inline Halfedge BoundaryLoopAdjacentHalfedgeNavigator::getCurrent() const { return currE; }

inline void BoundaryLoopAdjacentEdgeNavigator::advance() { currE = currE.next(); }
inline bool BoundaryLoopAdjacentEdgeNavigator::isValid() const { return true; }
inline Edge BoundaryLoopAdjacentEdgeNavigator::getCurrent() const { return currE.edge(); }

} // namespace halfedge_mesh
} // namespace geometrycentral
