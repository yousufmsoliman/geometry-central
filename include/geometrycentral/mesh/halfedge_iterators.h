#pragma once

#include <iterator>

#include "geometrycentral/mesh/halfedge_mesh.h"

namespace geometrycentral {
namespace halfedge_mesh {


// NOTE: These iterators are not STL compliant (so you can't use them with <algorithm> and friends). This is mainly
// becuase the STL notion of iterators seems to strongly imply each "container" has exactly one set of data to be
// iterated over. Obviously we break this here, we don't even have containers.

// The excessive inlining throughout these iterators was chosen after some halfhearted performance testing. When first
// implemented, the initial functions seemed to be a factor of 4 slower than the equivalent do{} while() loops. Now,
// they're are 1.0x-1.5x the cost.

// TODO add const interators across the board ?


// ==========================================================
// ================     Base  Iterator     ==================
// ==========================================================


// == Base range iterator
// All navigation iterators have the form "advance through elements, stopping when we get back where we started". The
// two classes below encapsulate that functionality, allowing us to just specify the element type, advancer, and "valid"
// function for each. The remaining boilerplates are generated. The template argument should be a class <N> with the
// following functionality:
//
//    - type Etype: the type of elements returned by the iterator (eg Vertex)
//    - member Etype currE: the current element pointed to by the iterator
//    - advance(): advance the iterator once (updating currE). Note that if some should be skipped, isValid will handle.
//    - isValid(): reports true if Etype should be returned by the iterator (not needed for most: always True)
//
// The helper classes will construct an instance of N, call advance() until isValid(), return currE. When iterator++ is
// invoked, the wrapper will call advance() until isValid() is satisfied. iterator== is implemented by comparing
// N::currE, in addition to a justStarted member which is set to false after the first increment. Probably best
// understood by reading the examples below.
template <typename N>
class NavigationIteratorBase {

public:
  NavigationIteratorBase(HalfedgeMesh* mesh_, typename N::Etype, bool justStarted_);
  const NavigationIteratorBase& operator++();
  bool operator==(const NavigationIteratorBase& other) const;
  bool operator!=(const NavigationIteratorBase& other) const;
  typename N::Etype operator*() const;

private:
  HalfedgeMesh* mesh;
  N state;
  bool justStarted; // our iterators are generally best expressed as "do-while" loops, so this is useful to distinguish
                    // begin() from end()
};

template <typename N>
class NavigationSetBase {
public:
  NavigationSetBase(HalfedgeMesh* mesh_, typename N::Etype);
  NavigationIteratorBase<N> begin() const;
  NavigationIteratorBase<N> end() const;

private:
  HalfedgeMesh* mesh;
  size_t iStart, iEnd;
};


// ==========================================================
// ================    Vertex Iterators    ==================
// ==========================================================

// Adjacent incoming halfedges
struct VertexIncomingHalfedgeNavigator {
  void advance();
  bool isValid() const;
  typedef Halfedge Etype;
  Etype currE;
};
typedef NavigationSetBase<VertexIncomingHalfedgeNavigator> VertexIncomingHalfedgeSet;

class VertexIncomingHalfedgeSet;
class VertexOutgoingHalfedgeSet;
class VertexAdjacentVertexSet;
class VertexAdjacentFaceSet;
class VertexAdjacentEdgeSet;
class VertexAdjacentCornerSet;

// Iterate around all incoming halfedges (both on the interior and the boundary
// of the domain)
class VertexIncomingHalfedgeIterator {
public:
  VertexIncomingHalfedgeIterator(Halfedge startingEdge, bool justStarted);
  const VertexIncomingHalfedgeIterator& operator++();
  bool operator==(const VertexIncomingHalfedgeIterator& other) const;
  bool operator!=(const VertexIncomingHalfedgeIterator& other) const;
  Halfedge operator*() const;

private:
  Halfedge currHe;
  bool justStarted;
};
class VertexIncomingHalfedgeSet {
public:
  VertexIncomingHalfedgeSet(Halfedge he);
  VertexIncomingHalfedgeIterator begin();
  VertexIncomingHalfedgeIterator end();

private:
  Halfedge firstHe;
};

// Iterate around all outgoing halfedges (both on the interior and the boundary
// of the domain)
class VertexOutgoingHalfedgeIterator {
public:
  VertexOutgoingHalfedgeIterator(Halfedge startingEdge, bool justStarted);
  const VertexOutgoingHalfedgeIterator& operator++();
  bool operator==(const VertexOutgoingHalfedgeIterator& other) const;
  bool operator!=(const VertexOutgoingHalfedgeIterator& other) const;
  Halfedge operator*() const;
  // Halfedge operator-> () const;

private:
  Halfedge currHe;
  bool justStarted;
};
class VertexOutgoingHalfedgeSet {
public:
  VertexOutgoingHalfedgeSet(Halfedge he);
  VertexOutgoingHalfedgeIterator begin();
  VertexOutgoingHalfedgeIterator end();

private:
  Halfedge firstHe;
};


// Iterate around adjacent vertices
class VertexAdjacentVertexIterator {
public:
  VertexAdjacentVertexIterator(Halfedge startingEdge, bool justStarted);
  const VertexAdjacentVertexIterator& operator++();
  bool operator==(const VertexAdjacentVertexIterator& other) const;
  bool operator!=(const VertexAdjacentVertexIterator& other) const;
  Vertex operator*() const;
  // Vertex operator-> () const;

private:
  Halfedge currHe;
  bool justStarted;
};
class VertexAdjacentVertexSet {
public:
  VertexAdjacentVertexSet(Halfedge he);
  VertexAdjacentVertexIterator begin();
  VertexAdjacentVertexIterator end();

private:
  Halfedge firstHe;
};

// Iterate around adjacent (real) faces
class VertexAdjacentFaceIterator {
public:
  VertexAdjacentFaceIterator(Halfedge startingEdge, bool justStarted);
  const VertexAdjacentFaceIterator& operator++();
  bool operator==(const VertexAdjacentFaceIterator& other) const;
  bool operator!=(const VertexAdjacentFaceIterator& other) const;
  Face operator*() const;
  // Face operator-> () const;

private:
  Halfedge currHe;
  bool justStarted;
};
class VertexAdjacentFaceSet {
public:
  VertexAdjacentFaceSet(Halfedge he);
  VertexAdjacentFaceIterator begin();
  VertexAdjacentFaceIterator end();

private:
  Halfedge firstHe;
};

// Iterate around adjacent edges
class VertexAdjacentEdgeIterator {
public:
  VertexAdjacentEdgeIterator(Halfedge startingEdge, bool justStarted);
  const VertexAdjacentEdgeIterator& operator++();
  bool operator==(const VertexAdjacentEdgeIterator& other) const;
  bool operator!=(const VertexAdjacentEdgeIterator& other) const;
  Edge operator*() const;
  // Edge operator-> () const;

private:
  Halfedge currHe;
  bool justStarted;
};
class VertexAdjacentEdgeSet {
public:
  VertexAdjacentEdgeSet(Halfedge he);
  VertexAdjacentEdgeIterator begin();
  VertexAdjacentEdgeIterator end();

private:
  Halfedge firstHe;
};

// Iterate around all adjacent corners
class VertexAdjacentCornerIterator {
public:
  VertexAdjacentCornerIterator(Halfedge startingEdge, bool justStarted);
  const VertexAdjacentCornerIterator& operator++();
  bool operator==(const VertexAdjacentCornerIterator& other) const;
  bool operator!=(const VertexAdjacentCornerIterator& other) const;
  Corner operator*() const;
  // Halfedge operator-> () const;

private:
  Halfedge currHe;
  bool justStarted;
};
class VertexAdjacentCornerSet {
public:
  VertexAdjacentCornerSet(Halfedge he);
  VertexAdjacentCornerIterator begin();
  VertexAdjacentCornerIterator end();

private:
  Halfedge firstHe;
};

// ==========================================================
// ================     Face Iterators     ==================
// ==========================================================

// Iterate around adjacent halfedges
class FaceAdjacentHalfedgeIterator {
public:
  FaceAdjacentHalfedgeIterator(Halfedge startingEdge, bool justStarted);
  const FaceAdjacentHalfedgeIterator& operator++();
  bool operator==(const FaceAdjacentHalfedgeIterator& other) const;
  bool operator!=(const FaceAdjacentHalfedgeIterator& other) const;
  Halfedge operator*() const;

private:
  Halfedge currHe;
  bool justStarted;
};
class FaceAdjacentHalfedgeSet {
public:
  FaceAdjacentHalfedgeSet(Halfedge he);
  FaceAdjacentHalfedgeIterator begin();
  FaceAdjacentHalfedgeIterator end();

private:
  Halfedge firstHe;
};

// Iterate around adjacent vertices
class FaceAdjacentVertexIterator {
public:
  FaceAdjacentVertexIterator(Halfedge startingEdge, bool justStarted);
  const FaceAdjacentVertexIterator& operator++();
  bool operator==(const FaceAdjacentVertexIterator& other) const;
  bool operator!=(const FaceAdjacentVertexIterator& other) const;
  Vertex operator*() const;

private:
  Halfedge currHe;
  bool justStarted;
};
class FaceAdjacentVertexSet {
public:
  FaceAdjacentVertexSet(Halfedge he);
  FaceAdjacentVertexIterator begin();
  FaceAdjacentVertexIterator end();

private:
  Halfedge firstHe;
};

// Iterate around adjacent edges
class FaceAdjacentEdgeIterator {
public:
  FaceAdjacentEdgeIterator(Halfedge startingEdge, bool justStarted);
  const FaceAdjacentEdgeIterator& operator++();
  bool operator==(const FaceAdjacentEdgeIterator& other) const;
  bool operator!=(const FaceAdjacentEdgeIterator& other) const;
  Edge operator*() const;

private:
  Halfedge currHe;
  bool justStarted;
};
class FaceAdjacentEdgeSet {
public:
  FaceAdjacentEdgeSet(Halfedge he);
  FaceAdjacentEdgeIterator begin();
  FaceAdjacentEdgeIterator end();

private:
  Halfedge firstHe;
};

// Iterate around adjacent (real) faces
class FaceAdjacentFaceIterator {
public:
  FaceAdjacentFaceIterator(Halfedge startingEdge, bool justStarted);
  const FaceAdjacentFaceIterator& operator++();
  bool operator==(const FaceAdjacentFaceIterator& other) const;
  bool operator!=(const FaceAdjacentFaceIterator& other) const;
  Face operator*() const;

private:
  Halfedge currHe;
  bool justStarted;
};
class FaceAdjacentFaceSet {
public:
  FaceAdjacentFaceSet(Halfedge he);
  FaceAdjacentFaceIterator begin();
  FaceAdjacentFaceIterator end();

private:
  Halfedge firstHe;
};

// Iterate around all adjacent corners
class FaceAdjacentCornerIterator {
public:
  FaceAdjacentCornerIterator(Halfedge startingEdge, bool justStarted_);
  const FaceAdjacentCornerIterator& operator++();
  bool operator==(const FaceAdjacentCornerIterator& other) const;
  bool operator!=(const FaceAdjacentCornerIterator& other) const;
  Corner operator*() const;
  // Halfedge operator-> () const;

private:
  Halfedge currHe;
  bool justStarted;
};
class FaceAdjacentCornerSet {
public:
  FaceAdjacentCornerSet(Halfedge he);
  FaceAdjacentCornerIterator begin();
  FaceAdjacentCornerIterator end();

private:
  Halfedge firstHe;
};


// ==========================================================
// ================     Face Iterators     ==================
// ==========================================================

// Iterate around adjacent halfedges
class BoundaryLoopAdjacentHalfedgeIterator {
public:
  BoundaryLoopAdjacentHalfedgeIterator(Halfedge startingEdge, bool justStarted);
  const BoundaryLoopAdjacentHalfedgeIterator& operator++();
  bool operator==(const BoundaryLoopAdjacentHalfedgeIterator& other) const;
  bool operator!=(const BoundaryLoopAdjacentHalfedgeIterator& other) const;
  Halfedge operator*() const;

private:
  Halfedge currHe;
  bool justStarted;
};
class BoundaryLoopAdjacentHalfedgeSet {
public:
  BoundaryLoopAdjacentHalfedgeSet(Halfedge he);
  BoundaryLoopAdjacentHalfedgeIterator begin();
  BoundaryLoopAdjacentHalfedgeIterator end();

private:
  Halfedge firstHe;
};

// Iterate around adjacent vertices
class BoundaryLoopAdjacentVertexIterator {
public:
  BoundaryLoopAdjacentVertexIterator(Halfedge startingEdge, bool justStarted);
  const BoundaryLoopAdjacentVertexIterator& operator++();
  bool operator==(const BoundaryLoopAdjacentVertexIterator& other) const;
  bool operator!=(const BoundaryLoopAdjacentVertexIterator& other) const;
  Vertex operator*() const;

private:
  Halfedge currHe;
  bool justStarted;
};
class BoundaryLoopAdjacentVertexSet {
public:
  BoundaryLoopAdjacentVertexSet(Halfedge he);
  BoundaryLoopAdjacentVertexIterator begin();
  BoundaryLoopAdjacentVertexIterator end();

private:
  Halfedge firstHe;
};

// Iterate around adjacent edges
class BoundaryLoopAdjacentEdgeIterator {
public:
  BoundaryLoopAdjacentEdgeIterator(Halfedge startingEdge, bool justStarted);
  const BoundaryLoopAdjacentEdgeIterator& operator++();
  bool operator==(const BoundaryLoopAdjacentEdgeIterator& other) const;
  bool operator!=(const BoundaryLoopAdjacentEdgeIterator& other) const;
  Edge operator*() const;

private:
  Halfedge currHe;
  bool justStarted;
};
class BoundaryLoopAdjacentEdgeSet {
public:
  BoundaryLoopAdjacentEdgeSet(Halfedge he);
  BoundaryLoopAdjacentEdgeIterator begin();
  BoundaryLoopAdjacentEdgeIterator end();

private:
  Halfedge firstHe;
};


} // namespace halfedge_mesh
} // namespace geometrycentral
