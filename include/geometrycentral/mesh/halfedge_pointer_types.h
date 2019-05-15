#pragma once

#include <cstddef>
#include <functional>
#include <iostream>

namespace geometrycentral {

// === Types and inline methods for the halfedge mesh pointer and datatypes

// === Forward declare iterator set types (pointers may need to return these to
// support range-based for loops)
class VertexIncomingHalfedgeSet;
class VertexIncomingInteriorHalfedgeSet;
class VertexOutgoingHalfedgeSet;
class VertexOutgoingInteriorHalfedgeSet;
class VertexAdjacentVertexSet;
class VertexAdjacentFaceSet;
class VertexAdjacentEdgeSet;
class VertexAdjacentCornerSet;
class FaceAdjacentHalfedgeSet;
class FaceAdjacentVertexSet;
class FaceAdjacentEdgeSet;
class FaceAdjacentFaceSet;
class FaceAdjacentCornerSet;

// === Pointer types for mesh elements
class Halfedge;
class Corner;
class Vertex;
class Edge;
class Face;

// The dynamic variants are automatically updated when the mesh is compressed, (the standard variants are invalidated).
class DynamicHalfedge;
class DynamicVertex;
class DynamicEdge;
class DynamicFace;
class DynamicBoundaryLoop;

// Halfedge
class Halfedge {
public:
  Halfedge(); // defaults to null
  Halfedge(Halfedge* ptr);
  Halfedge(DynamicHalfedge he);

  // Connectivity
  Halfedge twin() const;
  Halfedge next() const;
  Halfedge prev() const;
  Vertex vertex() const;
  Edge edge() const;
  Face face() const;
  Corner corner() const;

  // Properties
  bool isReal() const;

  // Accessors
  Halfedge operator*();
  Halfedge operator*() const;
  Halfedge* operator->();
  const Halfedge* operator->() const;

  // Comparators
  bool operator==(const Halfedge& other) const;
  bool operator!=(const Halfedge& other) const;
  bool operator>(const Halfedge& other) const;
  bool operator>=(const Halfedge& other) const;
  bool operator<(const Halfedge& other) const;
  bool operator<=(const Halfedge& other) const;

  // The dynamic equivalent
  typedef DynamicHalfedge DynamicType;

protected:
  size_t ind = INVALID_IND;
  HalfedgeMesh* mesh;

  friend class HalfedgeMesh;
  friend class DynamicHalfedge;
  friend std::ostream& operator<<(std::ostream& output, const Halfedge& he);
  friend struct std::hash<Halfedge>;
};
std::ostream& operator<<(std::ostream& output, const Halfedge& he);

class DynamicHalfedge {
public:
  DynamicHalfedge(){};
  DynamicHalfedge(Halfedge* ptr, HalfedgeMesh* mesh);
  DynamicHalfedge(Halfedge ptr, HalfedgeMesh* mesh);

  DynamicHalfedge twin() const;
  DynamicHalfedge next() const;
  DynamicVertex vertex() const;
  DynamicEdge edge() const;
  DynamicFace face() const;

  // Despite the name, should not be used as a [0,N) index.
  size_t getInd() const;

  // Comparators
  bool operator==(const DynamicHalfedge& other) const;
  bool operator!=(const DynamicHalfedge& other) const;
  bool operator>(const DynamicHalfedge& other) const;
  bool operator>=(const DynamicHalfedge& other) const;
  bool operator<(const DynamicHalfedge& other) const;
  bool operator<=(const DynamicHalfedge& other) const;

  // The static equivalent
  typedef Halfedge StaticType;

private:
  friend class Halfedge;
  friend class HalfedgeMesh;
  friend struct std::hash<DynamicHalfedge>;

  size_t ind = -1;
  HalfedgeMesh* mesh = nullptr;
  // std::list<DynamicHalfedge*>::iterator listIt;
};

enum class HalfedgeSetType { Real, Imaginary, All };
class HalfedgeRangeIterator {
public:
  HalfedgeRangeIterator(Halfedge startingHalfedge, HalfedgeSetType type_, Halfedge end_);
  const HalfedgeRangeIterator& operator++();
  bool operator==(const HalfedgeRangeIterator& other) const;
  bool operator!=(const HalfedgeRangeIterator& other) const;
  Halfedge operator*() const;

private:
  Halfedge currHalfedge;
  HalfedgeSetType type;
  Halfedge end; // unfortunately needed to respect type option
};
class HalfedgeSet {
public:
  HalfedgeSet(Halfedge beginptr_, Halfedge endptr_, HalfedgeSetType type_);
  HalfedgeRangeIterator begin();
  HalfedgeRangeIterator end();

private:
  Halfedge beginptr, endptr;
  HalfedgeSetType type;
};

// Corner
class Corner {
public:
  Corner(); // defaults to nullptr
  Corner(Halfedge* ptr);

  // Connectivity
  Corner next() const;
  Corner prev() const;
  Halfedge halfedge() const;
  Vertex vertex() const;
  Face face() const;

  // Accessors
  Halfedge operator*();
  Halfedge operator*() const;
  Halfedge* operator->();
  const Halfedge* operator->() const;

  // Comparators
  bool operator==(const Corner& other) const;
  bool operator!=(const Corner& other) const;
  bool operator>(const Corner& other) const;
  bool operator>=(const Corner& other) const;
  bool operator<(const Corner& other) const;
  bool operator<=(const Corner& other) const;

  // Null comparators
  bool operator==(void* n) const;
  bool operator!=(void* n) const;
  bool operator>(void* n) const;
  bool operator>=(void* n) const;
  bool operator<(void* n) const;
  bool operator<=(void* n) const;

  // Arithmetic
  unsigned int operator-(const Corner& other) const;
  Corner& operator++();
  Corner operator++(int);
  Corner& operator--();
  Corner operator--(int);

  // The dynamic equivalent
  // TODO
  typedef DynamicHalfedge DynamicType;

protected:
  Halfedge* ptr = nullptr;

  friend class HalfedgeMesh;
};
class CornerRangeIterator {
public:
  CornerRangeIterator(Corner startingCorner, Corner end);
  const CornerRangeIterator& operator++();
  bool operator==(const CornerRangeIterator& other) const;
  bool operator!=(const CornerRangeIterator& other) const;
  Corner operator*() const;

private:
  Corner currCorner;
  Corner end;
};
class CornerSet {
public:
  CornerSet(Corner beginptr_, Corner endptr_);
  CornerRangeIterator begin();
  CornerRangeIterator end();

private:
  Corner beginptr, endptr;
};

// Vertex
class Vertex {
public:
  Vertex(); // defaults to nullptr
  Vertex(Vertex* ptr);
  Vertex(DynamicVertex v);

  // Connectivity
  Halfedge halfedge() const;
  Corner corner() const;


  // Properties
  bool isBoundary() const;
  unsigned int degree();

  // Iterators
  VertexIncomingHalfedgeSet incomingHalfedges();
  VertexOutgoingHalfedgeSet outgoingHalfedges();
  VertexIncomingInteriorHalfedgeSet incomingInteriorHalfedges();
  VertexOutgoingInteriorHalfedgeSet outgoingInteriorHalfedges();
  VertexAdjacentVertexSet adjacentVertices();
  VertexAdjacentFaceSet adjacentFaces();
  VertexAdjacentEdgeSet adjacentEdges();
  VertexAdjacentCornerSet adjacentCorners();

  // Accessors
  Vertex operator*();
  Vertex operator*() const;
  Vertex* operator->();
  const Vertex* operator->() const;

  // Comparators
  bool operator==(const Vertex& other) const;
  bool operator!=(const Vertex& other) const;
  bool operator>(const Vertex& other) const;
  bool operator>=(const Vertex& other) const;
  bool operator<(const Vertex& other) const;
  bool operator<=(const Vertex& other) const;

  // Null comparators
  bool operator==(void* n) const;
  bool operator!=(void* n) const;
  bool operator>(void* n) const;
  bool operator>=(void* n) const;
  bool operator<(void* n) const;
  bool operator<=(void* n) const;

  // Arithmetic
  unsigned int operator-(const Vertex& other) const;
  Vertex& operator++();
  Vertex operator++(int);
  Vertex& operator--();
  Vertex operator--(int);

  // The dynamic equivalent
  typedef DynamicVertex DynamicType;

protected:
  Vertex* ptr = nullptr;

  friend class HalfedgeMesh;
  friend class DynamicVertex;
  friend std::ostream& operator<<(std::ostream& output, const Vertex& v);
  friend struct std::hash<Vertex>;
};
std::ostream& operator<<(std::ostream& output, const Vertex& v);

class DynamicVertex {
public:
  DynamicVertex(){};
  DynamicVertex(Vertex* ptr, HalfedgeMesh* mesh);
  DynamicVertex(Vertex ptr, HalfedgeMesh* mesh);

  // Connectivity
  DynamicHalfedge halfedge() const;

  // Despite the name, should not be used as a [0,N) index.
  size_t getInd() const;

  // Comparators
  bool operator==(const DynamicVertex& other) const;
  bool operator!=(const DynamicVertex& other) const;
  bool operator>(const DynamicVertex& other) const;
  bool operator>=(const DynamicVertex& other) const;
  bool operator<(const DynamicVertex& other) const;
  bool operator<=(const DynamicVertex& other) const;

  // The static equivalent
  typedef Vertex StaticType;

private:
  friend class Vertex;
  friend class HalfedgeMesh;
  friend struct std::hash<DynamicVertex>;

  size_t ind = -1;
  HalfedgeMesh* mesh = nullptr;
  // std::list<DynamicVertex*>::iterator listIt;
};

class VertexRangeIterator {
public:
  VertexRangeIterator(Vertex startingVertex, Vertex end);
  const VertexRangeIterator& operator++();
  bool operator==(const VertexRangeIterator& other) const;
  bool operator!=(const VertexRangeIterator& other) const;
  Vertex operator*() const;

private:
  Vertex currVertex;
  Vertex end;
};
class VertexSet {
public:
  VertexSet(Vertex beginptr_, Vertex endptr_);
  VertexRangeIterator begin();
  VertexRangeIterator end();

private:
  Vertex beginptr, endptr;
};

// Edge
class Edge {
public:
  Edge(); // defaults to nullptr
  Edge(Edge* ptr);
  Edge(DynamicEdge e);

  // Connectivity
  Halfedge halfedge() const;

  // Properties
  bool isBoundary() const;

  // Accessors
  Edge operator*();
  Edge operator*() const;
  Edge* operator->();
  const Edge* operator->() const;

  // Comparators
  bool operator==(const Edge& other) const;
  bool operator!=(const Edge& other) const;
  bool operator>(const Edge& other) const;
  bool operator>=(const Edge& other) const;
  bool operator<(const Edge& other) const;
  bool operator<=(const Edge& other) const;

  // Null comparators
  bool operator==(void* n) const;
  bool operator!=(void* n) const;
  bool operator>(void* n) const;
  bool operator>=(void* n) const;
  bool operator<(void* n) const;
  bool operator<=(void* n) const;

  // Arithmetic
  unsigned int operator-(const Edge& other) const;
  Edge& operator++();
  Edge operator++(int);
  Edge& operator--();
  Edge operator--(int);

  // The dynamic equivalent
  typedef DynamicEdge DynamicType;

protected:
  Edge* ptr = nullptr;

  friend class HalfedgeMesh;
  friend class DynamicEdge;
  friend std::ostream& operator<<(std::ostream& output, const Edge& e);
  friend struct std::hash<Edge>;
};
std::ostream& operator<<(std::ostream& output, const Edge& e);

class DynamicEdge {
public:
  DynamicEdge(){};
  DynamicEdge(Edge* ptr, HalfedgeMesh* mesh);
  DynamicEdge(Edge ptr, HalfedgeMesh* mesh);

  // Connectivity
  DynamicHalfedge halfedge() const;

  // Despite the name, should not be used as a [0,N) index.
  size_t getInd() const;

  // Comparators
  bool operator==(const DynamicEdge& other) const;
  bool operator!=(const DynamicEdge& other) const;
  bool operator>(const DynamicEdge& other) const;
  bool operator>=(const DynamicEdge& other) const;
  bool operator<(const DynamicEdge& other) const;
  bool operator<=(const DynamicEdge& other) const;

  // The static equivalent
  typedef Edge StaticType;

private:
  friend class Edge;
  friend class HalfedgeMesh;
  friend struct std::hash<DynamicEdge>;

  size_t ind = -1;
  HalfedgeMesh* mesh = nullptr;
  // std::list<DynamicEdge*>::iterator listIt;
};

class EdgeRangeIterator {
public:
  EdgeRangeIterator(Edge startingEdge, Edge end);
  const EdgeRangeIterator& operator++();
  bool operator==(const EdgeRangeIterator& other) const;
  bool operator!=(const EdgeRangeIterator& other) const;
  Edge operator*() const;

private:
  Edge currEdge;
  Edge end;
};
class EdgeSet {
public:
  EdgeSet(Edge beginptr_, Edge endptr_);
  EdgeRangeIterator begin();
  EdgeRangeIterator end();

private:
  Edge beginptr, endptr;
};

// NOTE: Triangle is merely a helper class for the
// method Face::triangulation(); it is NOT the
// standard representation for faces of a HalfedgeMesh.
struct Triangle {
  Vertex vertex[3];
  Vertex& operator[](size_t index) { return vertex[index]; }
};

// Face
class Face {
public:
  Face(); // defaults to nullptr
  Face(Face* ptr);
  Face(DynamicFace f);

  // Connectivity
  Halfedge halfedge() const;
  Corner corner() const;


  // Utility
  std::vector<Triangle> triangulation();

  // Properties
  unsigned int degree();
  bool isBoundary() const;
  bool isReal() const;

  // Iterators
  FaceAdjacentHalfedgeSet adjacentHalfedges();
  FaceAdjacentVertexSet adjacentVertices();
  FaceAdjacentFaceSet adjacentFaces();
  FaceAdjacentEdgeSet adjacentEdges();
  FaceAdjacentCornerSet adjacentCorners();

  // Accessors
  Face operator*();
  Face operator*() const;
  Face* operator->();
  const Face* operator->() const;

  // Comparators
  bool operator==(const Face& other) const;
  bool operator!=(const Face& other) const;
  bool operator>(const Face& other) const;
  bool operator>=(const Face& other) const;
  bool operator<(const Face& other) const;
  bool operator<=(const Face& other) const;

  // Null comparators
  bool operator==(void* n) const;
  bool operator!=(void* n) const;
  bool operator>(void* n) const;
  bool operator>=(void* n) const;
  bool operator<(void* n) const;
  bool operator<=(void* n) const;

  // Arithmetic
  unsigned int operator-(const Face& other) const;
  Face& operator++();
  Face operator++(int);
  Face& operator--();
  Face operator--(int);

  // The dynamic equivalent
  typedef DynamicFace DynamicType;

protected:
  Face* ptr = nullptr;

  friend class HalfedgeMesh;
  friend class DynamicFace;
  friend std::ostream& operator<<(std::ostream& output, const Face& f);
  friend struct std::hash<Face>;
};
std::ostream& operator<<(std::ostream& output, const Face& f);

class DynamicFace {
public:
  DynamicFace(){};
  DynamicFace(Face* ptr, HalfedgeMesh* mesh);
  DynamicFace(Face ptr, HalfedgeMesh* mesh);

  // Connectivity
  DynamicHalfedge halfedge() const;

  // Despite the name, should not be used as a [0,N) index.
  size_t getInd() const;

  // Comparators
  bool operator==(const DynamicFace& other) const;
  bool operator!=(const DynamicFace& other) const;
  bool operator>(const DynamicFace& other) const;
  bool operator>=(const DynamicFace& other) const;
  bool operator<(const DynamicFace& other) const;
  bool operator<=(const DynamicFace& other) const;

  // The static equivalent
  typedef Face StaticType;

private:
  friend class Face;
  friend class HalfedgeMesh;
  friend struct std::hash<DynamicFace>;

  size_t ind = -1;
  HalfedgeMesh* mesh = nullptr;
  // std::list<DynamicFace*>::iterator listIt;
};

class FaceRangeIterator {
public:
  FaceRangeIterator(Face startingFace, Face end);
  const FaceRangeIterator& operator++();
  bool operator==(const FaceRangeIterator& other) const;
  bool operator!=(const FaceRangeIterator& other) const;
  Face operator*() const;

private:
  Face currFace;
  Face end;
};
class FaceSet {
public:
  FaceSet(Face beginptr_, Face endptr_);
  FaceRangeIterator begin();
  FaceRangeIterator end();

private:
  Face beginptr, endptr;
};

// Boundary (currently just a renaming of Face---if we wanted
// stronger type checking we could instead inherit from Face*)

typedef Face BoundaryLoop;
typedef Face BoundaryLoop;
typedef FaceSet BoundarySet;
typedef FaceRangeIterator BoundaryRangeIterator;


} // namespace geometrycentral
