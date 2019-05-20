#include "geometrycentral/mesh/halfedge_mesh.h"

// Implementation of various methods to support iterators in the halfedge mesh.
// Note that many functions are inlined, and just defined in
// halfedge_iterators.ipp

namespace geometrycentral {
namespace halfedge_mesh {

// ==========================================================
// ================     Base  Iterator     ==================
// ==========================================================

// Note: the class below may seem a bit weird, in that we call advance() lazily in operator* rather than proactively in
// operator++ (and occaisonally in operator++ in case operator* was never called). We do this to avoid repeatedly
// advancing to hte first valid element on std::end(set). The intuition is that the iterators will still be efficient
// even if advancing to the first valid element is expensive (though we do assume that calling isValid() is cheap).

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
inline typename N::Etype NavigationIteratorBase<N>::operator*() const {
  while (!state.isValid()) {
    state.advance();
  }
  return state.currE;
}

template <typename N>
NavigationSetBase<N>::NavigationSetBase(HalfedgeMesh* mesh_, size_t iStart_, size_t iEnd_)
    : mesh(mesh_), iStart(iStart_), iEnd(iEnd_) {}

template <typename N>
inline NavigationIteratorBase<N> NavigationSetBase<N>::begin() const {
  return NavigationIteratorBase<N>(mesh, iStart, iEnd);
}

template <typename N>
inline NavigationIteratorBase<N> NavigationSetBase<N>::end() const {
  return NavigationIteratorBase<N>(mesh, iEnd, iEnd);
}


// ==========================================================
// ================    Vertex Iterators    ==================
// ==========================================================

// === Incoming halfedges (with imaginary)
VertexIncomingHalfedgeSet::VertexIncomingHalfedgeSet(Halfedge he) : firstHe(he) {}

VertexIncomingHalfedgeIterator::VertexIncomingHalfedgeIterator(Halfedge startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

// === Incoming halfedges (without imaginary)
VertexIncomingInteriorHalfedgeSet::VertexIncomingInteriorHalfedgeSet(Halfedge he) : firstHe(he) {}

VertexIncomingInteriorHalfedgeIterator::VertexIncomingInteriorHalfedgeIterator(Halfedge startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {
  // If startingEdge is a boundary edge, the first state of this iterator would
  // be a halfedge which should not be returned,
  // so advance the iterator until we find valid one.
  // TODO initialize the mesh so that the starting edge satisfies this property
  // by default

  while (!currHe.isReal()) {
    currHe = currHe.next().twin();
  }
}

// === Outgoing halfedges (with imaginary)
VertexOutgoingHalfedgeSet::VertexOutgoingHalfedgeSet(Halfedge he) : firstHe(he) {}

VertexOutgoingHalfedgeIterator::VertexOutgoingHalfedgeIterator(Halfedge startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

// === Outgoing halfedges (without imaginary)
VertexOutgoingInteriorHalfedgeSet::VertexOutgoingInteriorHalfedgeSet(Halfedge he) : firstHe(he) {}

VertexOutgoingInteriorHalfedgeIterator::VertexOutgoingInteriorHalfedgeIterator(Halfedge startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {
  // If startingEdge is a boundary edge, the first state of this iterator would
  // be a halfedge which should not be returned,
  // so advance the iterator until we find valid one.
  // TODO initialize the mesh so that the starting edge satisfies this property
  // by default

  while (!currHe.isReal()) {
    currHe = currHe.twin().next();
  }
}

// === Adjacent vertices
VertexAdjacentVertexSet::VertexAdjacentVertexSet(Halfedge he) : firstHe(he) {}

VertexAdjacentVertexIterator::VertexAdjacentVertexIterator(Halfedge startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

// === Adjacent faces
VertexAdjacentFaceSet::VertexAdjacentFaceSet(Halfedge he) : firstHe(he) {}

VertexAdjacentFaceIterator::VertexAdjacentFaceIterator(Halfedge startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

// === Adjacent edges
VertexAdjacentEdgeSet::VertexAdjacentEdgeSet(Halfedge he) : firstHe(he) {}

VertexAdjacentEdgeIterator::VertexAdjacentEdgeIterator(Halfedge startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

// === Adjacent corners
VertexAdjacentCornerSet::VertexAdjacentCornerSet(Halfedge he) : firstHe(he) {
  // Set firstHe to a real halfedge
  if (!firstHe.isReal()) firstHe = firstHe.twin().next();
}

VertexAdjacentCornerIterator::VertexAdjacentCornerIterator(Halfedge startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

// ==========================================================
// ================     Face Iterators     ==================
// ==========================================================

// === Adjacent halfedges
FaceAdjacentHalfedgeSet::FaceAdjacentHalfedgeSet(Halfedge he) : firstHe(he) {}

FaceAdjacentHalfedgeIterator::FaceAdjacentHalfedgeIterator(Halfedge startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

// === Adjacent vertices
FaceAdjacentVertexSet::FaceAdjacentVertexSet(Halfedge he) : firstHe(he) {}

FaceAdjacentVertexIterator::FaceAdjacentVertexIterator(Halfedge startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

// === Adjacent edges
FaceAdjacentEdgeSet::FaceAdjacentEdgeSet(Halfedge he) : firstHe(he) {}

FaceAdjacentEdgeIterator::FaceAdjacentEdgeIterator(Halfedge startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

// === Adjacent faces
FaceAdjacentFaceSet::FaceAdjacentFaceSet(Halfedge he) : firstHe(he) {}

FaceAdjacentFaceIterator::FaceAdjacentFaceIterator(Halfedge startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {
  // face->halfedge could point to a halfedge on the boundary, whose
  // corresponding face should
  // not be returned by this iterator. As such, advance the iterator during
  // construction until
  // it points to a non-border halfedge (or until it loops around).
  // An alternative would be to require that face->halfedge points to a
  // non-boundary halfedge,
  // but this feels like an unnatural requirement (for instance, it disallows a
  // single lonely
  // triangle as a mesh element).
  // As an aside, in the "lonely triangle" case, this code is correct but does
  // an extra
  // circuit of the face.
  // TODO can we initialize the mesh such that this loop is not needed (i.e.,
  // such that it's satisfied by default?)
  if (!currHe.twin().isReal()) {
    do {
      currHe = currHe.next();
    } while (!currHe.twin().isReal() && currHe != startingEdge);
  }
}

// === Adjacent corners
FaceAdjacentCornerSet::FaceAdjacentCornerSet(Halfedge he) : firstHe(he) {}

FaceAdjacentCornerIterator::FaceAdjacentCornerIterator(Halfedge startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

} // namespace halfedge_mesh
} // namespace geometrycentral
