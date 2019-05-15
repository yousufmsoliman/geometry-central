#pragma once

// Implementations for halfedge_mesh_types.ipp

namespace geometrycentral {

// Halfedge
inline Halfedge::Halfedge() : ptr(nullptr) {}
inline Halfedge::Halfedge(Halfedge* ptr_) : ptr(ptr_) {}
inline Halfedge::Halfedge(DynamicHalfedge ptr_) : Halfedge(ptr_.mesh->halfedge(ptr_.ind)) {}
inline Halfedge Halfedge::twin() const { return ptr->twin; }
inline Halfedge Halfedge::next() const { return ptr->next; }
inline Halfedge Halfedge::prev() const {
  Halfedge* h = ptr;
  while (h->next != ptr) {
    h = h->next;
  }
  return Halfedge{h};
}
inline Vertex Halfedge::vertex() const { return ptr->vertex; }
inline Edge Halfedge::edge() const { return ptr->edge; }
inline Face Halfedge::face() const { return ptr->face; }
inline Corner Halfedge::corner() const { return Corner{ptr}; }
inline bool Halfedge::isReal() const { return ptr->isReal; }
inline Halfedge Halfedge::operator*() { return *ptr; }
inline Halfedge Halfedge::operator*() const { return *ptr; }
inline Halfedge* Halfedge::operator->() { return ptr; }
inline const Halfedge* Halfedge::operator->() const { return ptr; }
inline bool Halfedge::operator==(const Halfedge& other) const { return ptr == other.ptr; }
inline bool Halfedge::operator!=(const Halfedge& other) const { return !(*this == other); }
inline bool Halfedge::operator>(const Halfedge& other) const { return ptr > other.ptr; }
inline bool Halfedge::operator>=(const Halfedge& other) const { return ptr >= other.ptr; }
inline bool Halfedge::operator<(const Halfedge& other) const { return ptr < other.ptr; }
inline bool Halfedge::operator<=(const Halfedge& other) const { return ptr <= other.ptr; }

inline ::std::ostream& operator<<(::std::ostream& output, const Halfedge& he) {
  output << "he_" << he.ptr;
  return output;
}

// FIXME dynamic pointers are currently broken -- need to support permutation callback for on compress()
inline DynamicHalfedge::DynamicHalfedge(Halfedge* ptr_, HalfedgeMesh* mesh_)
    : ind(mesh_->indexOf(ptr_)), mesh(mesh_) {}
inline DynamicHalfedge::DynamicHalfedge(Halfedge ptr_, HalfedgeMesh* mesh_)
    : ind(mesh_->indexOf(ptr_.ptr)), mesh(mesh_) {}
inline DynamicHalfedge DynamicHalfedge::twin() const {
  return DynamicHalfedge(Halfedge(*this).twin(), mesh);
}
inline DynamicHalfedge DynamicHalfedge::next() const {
  return DynamicHalfedge(Halfedge(*this).next(), mesh);
}
inline DynamicVertex DynamicHalfedge::vertex() const {
  return DynamicVertex(Halfedge(*this).vertex(), mesh);
}
inline DynamicEdge DynamicHalfedge::edge() const { return DynamicEdge(Halfedge(*this).edge(), mesh); }
inline DynamicFace DynamicHalfedge::face() const { return DynamicFace(Halfedge(*this).face(), mesh); }

inline bool DynamicHalfedge::operator==(const DynamicHalfedge& other) const { return ind == other.ind; }
inline bool DynamicHalfedge::operator!=(const DynamicHalfedge& other) const { return !(*this == other); }
inline bool DynamicHalfedge::operator>(const DynamicHalfedge& other) const { return ind > other.ind; }
inline bool DynamicHalfedge::operator>=(const DynamicHalfedge& other) const { return ind >= other.ind; }
inline bool DynamicHalfedge::operator<(const DynamicHalfedge& other) const { return ind < other.ind; }
inline bool DynamicHalfedge::operator<=(const DynamicHalfedge& other) const { return ind <= other.ind; }
inline size_t DynamicHalfedge::getInd() const { return ind; }

inline HalfedgeRangeIterator::HalfedgeRangeIterator(Halfedge startingHalfedge, HalfedgeSetType type_,
                                                          Halfedge end_)
    : currHalfedge(startingHalfedge), type(type_), end(end_) {
  // Advance to satisfy the type
  if (currHalfedge != end &&
      ((type == HalfedgeSetType::Real && !currHalfedge.isReal()) ||
       (type == HalfedgeSetType::Imaginary && currHalfedge.isReal()) || currHalfedge->isDead())) {
    this->operator++();
  }
}
inline const HalfedgeRangeIterator& HalfedgeRangeIterator::operator++() {
  currHalfedge++;

  // = Respect the 'type' option
  // Note that we always need to return if we fall off the end of the list, and must not dereference the pointer in that
  // case

  // All halfedges
  if (currHalfedge == end) {
    return *this;
  }

  if (type == HalfedgeSetType::All) {
    while (currHalfedge != end && currHalfedge->isDead()) {
      currHalfedge++;
    }
    return *this;
  }
  // Real only
  else if (type == HalfedgeSetType::Real) {
    while (currHalfedge != end && (!currHalfedge.isReal() || currHalfedge->isDead())) {
      currHalfedge++;
    }
    return *this;
  }
  // Imaginary only
  else /* imag */ {
    while (currHalfedge != end && (currHalfedge.isReal() || currHalfedge->isDead())) {
      currHalfedge++;
    }
    return *this;
  }
}
inline bool HalfedgeRangeIterator::operator==(const HalfedgeRangeIterator& other) const {
  return currHalfedge == other.currHalfedge;
}
inline bool HalfedgeRangeIterator::operator!=(const HalfedgeRangeIterator& other) const {
  return currHalfedge != other.currHalfedge;
}
inline Halfedge HalfedgeRangeIterator::operator*() const { return currHalfedge; }

inline HalfedgeSet::HalfedgeSet(Halfedge beginptr_, Halfedge endptr_, HalfedgeSetType type_)
    : beginptr(beginptr_), endptr(endptr_), type(type_) {}
inline HalfedgeRangeIterator HalfedgeSet::begin() { return HalfedgeRangeIterator(beginptr, type, endptr); }
inline HalfedgeRangeIterator HalfedgeSet::end() { return HalfedgeRangeIterator(endptr, type, endptr); }


// Corner
inline Corner::Corner() : ptr(nullptr) {}
inline Corner::Corner(Halfedge* ptr_) : ptr(ptr_) {}
inline Corner Corner::next() const { return halfedge().next().corner(); }
inline Corner Corner::prev() const { return halfedge().prev().corner(); }
inline Halfedge Corner::halfedge() const { return Halfedge{ptr}; }
inline Vertex Corner::vertex() const { return halfedge().prev().vertex(); }
inline Face Corner::face() const { return halfedge().face(); }
inline Halfedge Corner::operator*() { return *ptr; }
inline Halfedge Corner::operator*() const { return *ptr; }
inline Halfedge* Corner::operator->() { return ptr; }
inline const Halfedge* Corner::operator->() const { return ptr; }
inline bool Corner::operator==(const Corner& other) const { return ptr == other.ptr; }
inline bool Corner::operator!=(const Corner& other) const { return !(*this == other); }
inline bool Corner::operator>(const Corner& other) const { return ptr > other.ptr; }
inline bool Corner::operator>=(const Corner& other) const { return ptr >= other.ptr; }
inline bool Corner::operator<(const Corner& other) const { return ptr < other.ptr; }
inline bool Corner::operator<=(const Corner& other) const { return ptr <= other.ptr; }
inline bool Corner::operator==(void* n) const { return ptr == n; }
inline bool Corner::operator!=(void* n) const { return ptr != n; }
inline bool Corner::operator>(void* n) const { return ptr > n; }
inline bool Corner::operator>=(void* n) const { return ptr >= n; }
inline bool Corner::operator<(void* n) const { return ptr < n; }
inline bool Corner::operator<=(void* n) const { return ptr <= n; }
inline unsigned int Corner::operator-(const Corner& other) const {
#ifndef NDEBUG
  assert(ptr >= other.ptr && "Pointer subtraction must yield a nonnegative result");
#endif
  return ptr - other.ptr;
}
inline Corner& Corner::operator++() {
  ptr++;
  return *this;
}
inline Corner Corner::operator++(int) {
  ptr++;
  return Corner(ptr - 1);
}
inline Corner& Corner::operator--() {
  ptr--;
  return *this;
}
inline Corner Corner::operator--(int) {
  ptr--;
  return Corner(ptr + 1);
}

inline CornerRangeIterator::CornerRangeIterator(Corner startingCorner, Corner end_)
    : currCorner(startingCorner), end(end_) {
  if (currCorner != end && currCorner->isDead()) {
    this->operator++();
  }
}
inline const CornerRangeIterator& CornerRangeIterator::operator++() {
  currCorner++;
  while (currCorner != end && currCorner->isDead()) {
    currCorner++;
  }
  return *this;
}
inline bool CornerRangeIterator::operator==(const CornerRangeIterator& other) const {
  return currCorner == other.currCorner;
}
inline bool CornerRangeIterator::operator!=(const CornerRangeIterator& other) const {
  return currCorner != other.currCorner;
}
inline Corner CornerRangeIterator::operator*() const { return currCorner; }

inline CornerSet::CornerSet(Corner beginptr_, Corner endptr_) : beginptr(beginptr_), endptr(endptr_) {}
inline CornerRangeIterator CornerSet::begin() { return CornerRangeIterator(beginptr, endptr); }
inline CornerRangeIterator CornerSet::end() { return CornerRangeIterator(endptr, endptr); }

// Vertex
inline Vertex::Vertex() : ptr(nullptr) {}
inline Vertex::Vertex(Vertex* ptr_) : ptr(ptr_) {}
inline Vertex::Vertex(DynamicVertex ptr_) : Vertex(ptr_.mesh->vertex(ptr_.ind)) {}
inline Halfedge Vertex::halfedge() const { return ptr->halfedge; }
inline Corner Vertex::corner() const {
  Halfedge h = halfedge();
  if (!h.isReal()) h = h.twin().next();
  return h.next().corner();
}
inline bool Vertex::isBoundary() const { return ptr->isBoundary; }
inline VertexIncomingHalfedgeSet Vertex::incomingHalfedges() { return VertexIncomingHalfedgeSet(halfedge().twin()); }
inline VertexOutgoingHalfedgeSet Vertex::outgoingHalfedges() { return VertexOutgoingHalfedgeSet(halfedge()); }
inline VertexIncomingInteriorHalfedgeSet Vertex::incomingInteriorHalfedges() {
  return VertexIncomingInteriorHalfedgeSet(halfedge().twin());
}
inline VertexOutgoingInteriorHalfedgeSet Vertex::outgoingInteriorHalfedges() {
  return VertexOutgoingInteriorHalfedgeSet(halfedge());
}
inline VertexAdjacentVertexSet Vertex::adjacentVertices() { return VertexAdjacentVertexSet(halfedge().twin()); }
inline VertexAdjacentFaceSet Vertex::adjacentFaces() { return VertexAdjacentFaceSet(halfedge()); }
inline VertexAdjacentEdgeSet Vertex::adjacentEdges() { return VertexAdjacentEdgeSet(halfedge()); }
inline VertexAdjacentCornerSet Vertex::adjacentCorners() {
  Halfedge h = halfedge();
  if (!h.isReal()) h = h.twin().next();
  return VertexAdjacentCornerSet(h);
}
inline Vertex Vertex::operator*() { return *ptr; }
inline Vertex Vertex::operator*() const { return *ptr; }
inline Vertex* Vertex::operator->() { return ptr; }
inline const Vertex* Vertex::operator->() const { return ptr; }
inline bool Vertex::operator==(const Vertex& other) const { return ptr == other.ptr; }
inline bool Vertex::operator!=(const Vertex& other) const { return !(*this == other); }
inline bool Vertex::operator>(const Vertex& other) const { return ptr > other.ptr; }
inline bool Vertex::operator>=(const Vertex& other) const { return ptr >= other.ptr; }
inline bool Vertex::operator<(const Vertex& other) const { return ptr < other.ptr; }
inline bool Vertex::operator<=(const Vertex& other) const { return ptr <= other.ptr; }
inline bool Vertex::operator==(void* n) const { return ptr == n; }
inline bool Vertex::operator!=(void* n) const { return ptr != n; }
inline bool Vertex::operator>(void* n) const { return ptr > n; }
inline bool Vertex::operator>=(void* n) const { return ptr >= n; }
inline bool Vertex::operator<(void* n) const { return ptr < n; }
inline bool Vertex::operator<=(void* n) const { return ptr <= n; }
inline unsigned int Vertex::operator-(const Vertex& other) const {
#ifndef NDEBUG
  assert(ptr >= other.ptr && "Pointer subtraction must yield a nonnegative result");
#endif
  return ptr - other.ptr;
}
inline Vertex& Vertex::operator++() {
  ptr++;
  return *this;
}
inline Vertex Vertex::operator++(int) {
  ptr++;
  return Vertex(ptr - 1);
}
inline Vertex& Vertex::operator--() {
  ptr--;
  return *this;
}
inline Vertex Vertex::operator--(int) {
  ptr--;
  return Vertex(ptr + 1);
}

inline ::std::ostream& operator<<(::std::ostream& output, const Vertex& v) {
  output << "v_" << v.ptr;
  return output;
}

inline DynamicVertex::DynamicVertex(Vertex* ptr_, HalfedgeMesh* mesh_) : ind(mesh_->indexOf(ptr_)), mesh(mesh_) {}
inline DynamicVertex::DynamicVertex(Vertex ptr_, HalfedgeMesh* mesh_)
    : ind(mesh_->indexOf(ptr_.ptr)), mesh(mesh_) {}
inline DynamicHalfedge DynamicVertex::halfedge() const {
  return DynamicHalfedge(Vertex(*this).halfedge(), mesh);
}
inline bool DynamicVertex::operator==(const DynamicVertex& other) const { return ind == other.ind; }
inline bool DynamicVertex::operator!=(const DynamicVertex& other) const { return !(*this == other); }
inline bool DynamicVertex::operator>(const DynamicVertex& other) const { return ind > other.ind; }
inline bool DynamicVertex::operator>=(const DynamicVertex& other) const { return ind >= other.ind; }
inline bool DynamicVertex::operator<(const DynamicVertex& other) const { return ind < other.ind; }
inline bool DynamicVertex::operator<=(const DynamicVertex& other) const { return ind <= other.ind; }

inline size_t DynamicVertex::getInd() const { return ind; }

inline VertexRangeIterator::VertexRangeIterator(Vertex startingVertex, Vertex end_)
    : currVertex(startingVertex), end(end_) {
  if (currVertex != end && currVertex->isDead()) {
    this->operator++();
  }
}
inline const VertexRangeIterator& VertexRangeIterator::operator++() {
  currVertex++;
  while (currVertex != end && currVertex->isDead()) {
    currVertex++;
  }
  return *this;
}
inline bool VertexRangeIterator::operator==(const VertexRangeIterator& other) const {
  return currVertex == other.currVertex;
}
inline bool VertexRangeIterator::operator!=(const VertexRangeIterator& other) const {
  return currVertex != other.currVertex;
}
inline Vertex VertexRangeIterator::operator*() const { return currVertex; }

inline VertexSet::VertexSet(Vertex beginptr_, Vertex endptr_) : beginptr(beginptr_), endptr(endptr_) {}
inline VertexRangeIterator VertexSet::begin() { return VertexRangeIterator(beginptr, endptr); }
inline VertexRangeIterator VertexSet::end() { return VertexRangeIterator(endptr, endptr); }

// Edge
inline Edge::Edge() : ptr(nullptr) {}
inline Edge::Edge(Edge* ptr_) : ptr(ptr_) {}
inline Edge::Edge(DynamicEdge ptr_) : Edge(ptr_.mesh->edge(ptr_.ind)) {}
inline Halfedge Edge::halfedge() const { return ptr->halfedge; }
inline bool Edge::isBoundary() const { return ptr->isBoundary; }
inline Edge Edge::operator*() { return *ptr; }
inline Edge Edge::operator*() const { return *ptr; }
inline Edge* Edge::operator->() { return ptr; }
inline const Edge* Edge::operator->() const { return ptr; }
inline bool Edge::operator==(const Edge& other) const { return ptr == other.ptr; }
inline bool Edge::operator!=(const Edge& other) const { return !(*this == other); }
inline bool Edge::operator>(const Edge& other) const { return ptr > other.ptr; }
inline bool Edge::operator>=(const Edge& other) const { return ptr >= other.ptr; }
inline bool Edge::operator<(const Edge& other) const { return ptr < other.ptr; }
inline bool Edge::operator<=(const Edge& other) const { return ptr <= other.ptr; }
inline bool Edge::operator==(void* n) const { return ptr == n; }
inline bool Edge::operator!=(void* n) const { return ptr != n; }
inline bool Edge::operator>(void* n) const { return ptr > n; }
inline bool Edge::operator>=(void* n) const { return ptr >= n; }
inline bool Edge::operator<(void* n) const { return ptr < n; }
inline bool Edge::operator<=(void* n) const { return ptr <= n; }
inline unsigned int Edge::operator-(const Edge& other) const {
#ifndef NDEBUG
  assert(ptr >= other.ptr && "Pointer subtraction must yield a nonnegative result");
#endif
  return ptr - other.ptr;
}
inline Edge& Edge::operator++() {
  ptr++;
  return *this;
}
inline Edge Edge::operator++(int) {
  ptr++;
  return Edge(ptr - 1);
}
inline Edge& Edge::operator--() {
  ptr--;
  return *this;
}
inline Edge Edge::operator--(int) {
  ptr--;
  return Edge(ptr + 1);
}

inline ::std::ostream& operator<<(::std::ostream& output, const Edge& e) {
  output << "e_" << e.ptr;
  return output;
}

inline DynamicEdge::DynamicEdge(Edge* ptr_, HalfedgeMesh* mesh_) : ind(mesh_->indexOf(ptr_)), mesh(mesh_) {}
inline DynamicEdge::DynamicEdge(Edge ptr_, HalfedgeMesh* mesh_) : ind(mesh_->indexOf(ptr_.ptr)), mesh(mesh_) {}
inline DynamicHalfedge DynamicEdge::halfedge() const {
  return DynamicHalfedge(Edge(*this).halfedge(), mesh);
}
inline bool DynamicEdge::operator==(const DynamicEdge& other) const { return ind == other.ind; }
inline bool DynamicEdge::operator!=(const DynamicEdge& other) const { return !(*this == other); }
inline bool DynamicEdge::operator>(const DynamicEdge& other) const { return ind > other.ind; }
inline bool DynamicEdge::operator>=(const DynamicEdge& other) const { return ind >= other.ind; }
inline bool DynamicEdge::operator<(const DynamicEdge& other) const { return ind < other.ind; }
inline bool DynamicEdge::operator<=(const DynamicEdge& other) const { return ind <= other.ind; }
inline size_t DynamicEdge::getInd() const { return ind; }

inline EdgeRangeIterator::EdgeRangeIterator(Edge startingEdge, Edge end_)
    : currEdge(startingEdge), end(end_) {}
inline const EdgeRangeIterator& EdgeRangeIterator::operator++() {
  currEdge++;
  while (currEdge != end && currEdge->isDead()) {
    currEdge++;
  }
  return *this;
}
inline bool EdgeRangeIterator::operator==(const EdgeRangeIterator& other) const {
  return currEdge == other.currEdge;
}
inline bool EdgeRangeIterator::operator!=(const EdgeRangeIterator& other) const {
  return currEdge != other.currEdge;
}
inline Edge EdgeRangeIterator::operator*() const { return currEdge; }

inline EdgeSet::EdgeSet(Edge beginptr_, Edge endptr_) : beginptr(beginptr_), endptr(endptr_) {}
inline EdgeRangeIterator EdgeSet::begin() { return EdgeRangeIterator(beginptr, endptr); }
inline EdgeRangeIterator EdgeSet::end() { return EdgeRangeIterator(endptr, endptr); }

// Face
inline Face::Face() : ptr(nullptr) {}
inline Face::Face(Face* ptr_) : ptr(ptr_) {}
inline Face::Face(DynamicFace ptr_) : Face(ptr_.mesh->face(ptr_.ind)) {}
inline Halfedge Face::halfedge() const { return ptr->halfedge; }
inline Corner Face::corner() const { return halfedge().next().corner(); }
inline unsigned int Face::degree() {
  unsigned int k = 0;
  for (Halfedge h : adjacentHalfedges()) {
    k++;
  }
  return k;
}
inline bool Face::isBoundary() const { return ptr->isBoundary; }
inline bool Face::isReal() const { return ptr->isReal; }
inline FaceAdjacentHalfedgeSet Face::adjacentHalfedges() { return FaceAdjacentHalfedgeSet(halfedge()); }
inline FaceAdjacentVertexSet Face::adjacentVertices() { return FaceAdjacentVertexSet(halfedge()); }
inline FaceAdjacentFaceSet Face::adjacentFaces() { return FaceAdjacentFaceSet(halfedge()); }
inline FaceAdjacentEdgeSet Face::adjacentEdges() { return FaceAdjacentEdgeSet(halfedge()); }
inline FaceAdjacentCornerSet Face::adjacentCorners() { return FaceAdjacentCornerSet(halfedge()); }
inline Face Face::operator*() { return *ptr; }
inline Face Face::operator*() const { return *ptr; }
inline Face* Face::operator->() { return ptr; }
inline const Face* Face::operator->() const { return ptr; }
inline bool Face::operator==(const Face& other) const { return ptr == other.ptr; }
inline bool Face::operator!=(const Face& other) const { return !(*this == other); }
inline bool Face::operator>(const Face& other) const { return ptr > other.ptr; }
inline bool Face::operator>=(const Face& other) const { return ptr >= other.ptr; }
inline bool Face::operator<(const Face& other) const { return ptr < other.ptr; }
inline bool Face::operator<=(const Face& other) const { return ptr <= other.ptr; }
inline bool Face::operator==(void* n) const { return ptr == n; }
inline bool Face::operator!=(void* n) const { return ptr != n; }
inline bool Face::operator>(void* n) const { return ptr > n; }
inline bool Face::operator>=(void* n) const { return ptr >= n; }
inline bool Face::operator<(void* n) const { return ptr < n; }
inline bool Face::operator<=(void* n) const { return ptr <= n; }
inline unsigned int Face::operator-(const Face& other) const {
#ifndef NDEBUG
  assert(ptr >= other.ptr && "Pointer subtraction must yield a nonnegative result");
#endif
  return ptr - other.ptr;
}
inline Face& Face::operator++() {
  ptr++;
  return *this;
}
inline Face Face::operator++(int) {
  ptr++;
  return Face(ptr - 1);
}
inline Face& Face::operator--() {
  ptr--;
  return *this;
}
inline Face Face::operator--(int) {
  ptr--;
  return Face(ptr + 1);
}

inline ::std::ostream& operator<<(::std::ostream& output, const Face& f) {
  output << "f_" << f.ptr;
  return output;
}

inline DynamicFace::DynamicFace(Face* ptr_, HalfedgeMesh* mesh_) : ind(mesh_->indexOf(ptr_)), mesh(mesh_) {}
inline DynamicFace::DynamicFace(Face ptr_, HalfedgeMesh* mesh_) : ind(mesh_->indexOf(ptr_.ptr)), mesh(mesh_) {}
inline DynamicHalfedge DynamicFace::halfedge() const {
  return DynamicHalfedge(Face(*this).halfedge(), mesh);
}
inline bool DynamicFace::operator==(const DynamicFace& other) const { return ind == other.ind; }
inline bool DynamicFace::operator!=(const DynamicFace& other) const { return !(*this == other); }
inline bool DynamicFace::operator>(const DynamicFace& other) const { return ind > other.ind; }
inline bool DynamicFace::operator>=(const DynamicFace& other) const { return ind >= other.ind; }
inline bool DynamicFace::operator<(const DynamicFace& other) const { return ind < other.ind; }
inline bool DynamicFace::operator<=(const DynamicFace& other) const { return ind <= other.ind; }
inline size_t DynamicFace::getInd() const { return ind; }

inline FaceRangeIterator::FaceRangeIterator(Face startingFace, Face end_)
    : currFace(startingFace), end(end_) {
  if (currFace != end && currFace->isDead()) {
    this->operator++();
  }
}
inline const FaceRangeIterator& FaceRangeIterator::operator++() {
  currFace++;
  while (currFace != end && currFace->isDead()) {
    currFace++;
  }
  return *this;
}
inline bool FaceRangeIterator::operator==(const FaceRangeIterator& other) const {
  return currFace == other.currFace;
}
inline bool FaceRangeIterator::operator!=(const FaceRangeIterator& other) const {
  return currFace != other.currFace;
}
inline Face FaceRangeIterator::operator*() const { return currFace; }

inline FaceSet::FaceSet(Face beginptr_, Face endptr_) : beginptr(beginptr_), endptr(endptr_) {}
inline FaceRangeIterator FaceSet::begin() { return FaceRangeIterator(beginptr, endptr); }
inline FaceRangeIterator FaceSet::end() { return FaceRangeIterator(endptr, endptr); }

} // namespace geometrycentral

namespace std {

// == Hash standard pointers

template <>
struct hash<geometrycentral::Halfedge> {
  std::size_t operator()(const geometrycentral::Halfedge& he) const { return std::hash<size_t>{}(he.ptr->ID); }
};

template <>
struct hash<geometrycentral::Vertex> {
  std::size_t operator()(const geometrycentral::Vertex& v) const { return std::hash<size_t>{}(v.ptr->ID); }
};

template <>
struct hash<geometrycentral::Edge> {
  std::size_t operator()(const geometrycentral::Edge& e) const { return std::hash<size_t>{}(e.ptr->ID); }
};

template <>
struct hash<geometrycentral::Face> {
  std::size_t operator()(const geometrycentral::Face& f) const { return std::hash<size_t>{}(f.ptr->ID); }
};

// == Hash dynamic pointers

template <>
struct hash<geometrycentral::DynamicHalfedge> {
  std::size_t operator()(const geometrycentral::DynamicHalfedge& he) const { return std::hash<size_t>{}(he.ind); }
};

template <>
struct hash<geometrycentral::DynamicVertex> {
  std::size_t operator()(const geometrycentral::DynamicVertex& v) const { return std::hash<size_t>{}(v.ind); }
};

template <>
struct hash<geometrycentral::DynamicEdge> {
  std::size_t operator()(const geometrycentral::DynamicEdge& e) const { return std::hash<size_t>{}(e.ind); }
};

template <>
struct hash<geometrycentral::DynamicFace> {
  std::size_t operator()(const geometrycentral::DynamicFace& f) const { return std::hash<size_t>{}(f.ind); }
};
} // namespace std
