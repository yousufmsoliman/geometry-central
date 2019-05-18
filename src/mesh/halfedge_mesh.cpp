#include "geometrycentral/mesh/halfedge_mesh.h"

#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include <geometrycentral/utilities/disjoint_sets.h>
#include <geometrycentral/utilities/timing.h>

using std::cout;
using std::endl;

namespace geometrycentral {
namespace halfedge_mesh {


HalfedgeMesh::HalfedgeMesh() {}

// Helpers for below
namespace {

// Find an element in a sorted list
size_t halfedgeLookup(const std::vector<size_t>& compressedList, size_t target, size_t start, size_t end) {
  // Linear search is fast for small searches
  if (end - start < 20) {
    for (size_t i = start; i < end; i++) {
      if (compressedList[i] == target) {
        return i;
      }
    }
    return std::numeric_limits<size_t>::max();
  }
  // ...but we don't want to degrade to O(N^2) for really high valence vertices,
  // so fall back to a binary search
  else {
    auto loc = std::lower_bound(compressedList.begin() + start, compressedList.begin() + end, target);

    if (loc != (compressedList.begin() + end) && (target == *loc)) {
      return loc - compressedList.begin();
    } else {
      return std::numeric_limits<size_t>::max();
    }
  }
}

} // namespace

HalfedgeMesh::HalfedgeMesh(const std::vector<std::vector<size_t>>& polygons, bool verbose) {

  // Assumes that the input index set is dense. This sometimes isn't true of (eg) obj files floating around the
  // internet, so consider removing unused vertices first when reading from foreign sources.

  // Flatten in the input list and measure some element counts
  nFacesCount = polygons.size();
  nVerticesCount = 0;
  std::vector<size_t> flatFaces;
  std::vector<size_t> faceDegrees;
  faceDegrees.reserve(nFacesCount);
  for (auto poly : polygons) {
    GC_SAFETY_ASSERT(poly.size() >= 3, "faces must have degree >= 3");
    faceDegrees.push_back(poly.size());
    for (auto i : poly) {
      nVerticesCount = std::max(nVerticesCount, i);
      flatFaces.push_back(i);
    }
  }
  nVerticesCount++; // 0-based means count is max+1

  // Pre-allocate face and vertex arrays
  vHalfedge = std::vector<size_t>(nVerticesCount, INVALID_IND);
  fHalfedge = std::vector<size_t>(nVerticesCount, INVALID_IND);

  // Track halfedges which have already been created
  // TODO replace with compressed list for performance
  std::unordered_map<std::tuple<size_t, size_t>, size_t> createdHalfedges;
  auto createdHeLookup = [&](std::tuple<size_t, size_t> key) -> size_t& {
    if (createdHalfedges.find(key) == createdHalfedges.end()) {
      createdHalfedges[key] = INVALID_IND;
    }
    return createdHalfedges[key];
  };

  // Walk the faces, creating halfedges and hooking up pointers
  size_t iFlatHeStart = 0;
  for (size_t iFace = 0; iFace < nFaces; iFace++) {

    // Walk around this face
    size_t faceDegree = faceDegrees[iFace];
    size_t prevHeInd = INVALID_IND;
    size_t firstHeInd = INVALID_IND;
    for (size_t iFaceHe = 0; iFaceHe < faceDegree; iFaceHe++) {

      size_t indTail = iFace + iFaceHe;
      size_t indTip = iFace + (iFaceHe + 1) % faceDegree;

      // Get an index for this halfedge
      std::tuple<size_t, size_t> heKey{indTail, indTip};
      std::tuple<size_t, size_t> heTwinKey{indTip, indTail};
      size_t& halfedgeInd = createdHeLookup(heKey);

      // Some sanity checks
      GC_SAFETY_ASSERT(indTail != indTip,
                       "self-edge in face list " + std::to_string(indTail) + " -- " + std::to_string(indTip));
      GC_SAFETY_ASSERT(halfedgeInd == INVALID_IND,
                       "duplicate edge in list " + std::to_string(indTail) + " -- " + std::to_string(indTip));

      // Find the twin to check if the element is already created
      size_t twinInd = createdHeLookup(heTwinKey);
      if (twinInd == INVALID_IND) {
        // If we haven't seen the twin yet either, create a new edge
        // TODO use createHalfedge() or something here?
        halfedgeInd = nHalfedgesCount;
        nHalfedgesCount += 2;

        // Grow arrays to make space
        heNext.push_back(INVALID_IND);
        heNext.push_back(INVALID_IND);
        heVertex.push_back(INVALID_IND);
        heVertex.push_back(INVALID_IND);
        heFace.push_back(INVALID_IND);
        heFace.push_back(INVALID_IND);
      } else {
        // If the twin has already been created, we have an index for the halfedge
        halfedgeInd = heTwin(twinInd);
      }

      // Hook up a bunch of pointers
      heFace[halfedgeInd] = iFace;
      heVertex[halfedgeInd] = indTail;
      vHalfedge[indTail] = halfedgeInd;
      if (iFaceHe == 0) {
        fHalfedge[iFace] = halfedgeInd;
        firstHeInd = halfedgeInd;
      } else {
        heNext[prevHeInd] = halfedgeInd;
      }
      prevHeInd = halfedgeInd;
    }

    heNext[prevHeInd] = firstHeInd; // hook up the first next() pointer, which we missed in the loop above

    // Prepare to loop again
    iFlatHeStart += faceDegree;
  }


  // Ensure that each boundary neighborhood is either a disk or a half-disk. Harder to diagnose if we wait until the
  // boundary walk below.
#ifndef NGC_SAFTEY_CHECKS
  {
    std::vector<char> vertexOnBoundary(nVerticesCount, false);
    for (size_t iHe = 0; iHe < nHalfedgesCount; iHe++) {
      if (heNext[iHe] == INVALID_IND) {
        size_t v = heVertex[iHe];
        GC_SAFETY_ASSERT(!vertexOnBoundary[v],
                         "vertex " + std::to_string(iV) + " appears in more than one boundary loop");
        vertexOnBoundary[v] = true;
      }
    }
  }
#endif

  // == Resolve boundary loops
  nInteriorHalfedgesCount = nHalfedgesCount; // will decrement as we find exterior
  for (size_t iHe = 0; iHe < nHalfedgesCount; iHe++) {

    // If the face pointer is invalid, the halfedge must be along an unresolved boundary loop
    if (heFace[iHe] != INVALID_IND) continue;

    // Create the new boundary loop
    size_t boundaryLoopInd = nFacesCount + nBoundaryLoopsCount;
    fHalfedge.push_back(iHe);
    nBoundaryLoopsCount++;

    // = Walk around the loop (CW)
    size_t currHe = iHe;
    size_t prevHe = INVALID_IND;
    do {

      // The boundary loop is the face for these halfedges
      heFace[currHe] = boundaryLoopInd;

      // This isn't an interior halfedge.
      nInteriorHalfedgesCount--;

      // Advance to the next halfedge along the boundary
      prevHe = currHe;
      currHe = heTwin(heNext[heTwin(currHe)]);
      while (heFace[iHe] != INVALID_IND) {
        currHe = heNext[heTwin(currHe)];
      }

      // Set the next pointer around the boundary loop
      heNext[currHe] = prevHe;

    } while (currHe != iHe);
  }

  // SOMEDAY: could shrink_to_fit() std::vectors here, at the cost of a copy. What's preferable?

  // Set capacities and other properties
  nEdgesCount = nHalfedgesCount / 2;
  nVertexCapacityCount = nVerticesCount;
  nHalfedgeCapacityCount = nHalfedgesCount;
  nFaceCapacityCount = nFacesCount + nBoundaryLoopsCount;
  nVertexFillCount = nVerticesCount;
  nHalfedgeFillCount = nHalfedgesCount;
  nFaceFillCount = nFacesCapacity;       
  nBoundaryLoopFillCount = nBoundaryLoops;
  isCanonicalFlag = true;
  isCompressedFlag = true;
  

#ifndef NGC_SAFTEY_CHECKS
  { // Check that the input was manifold in the sense that each vertex has a single connected loop of faces around it.
    std::vector<char> halfedgeSeen(nHalfedgesCount, false);
    for (size_t iV = 0; iV < nVerticesCount; iV++) {

      // For each vertex, orbit around the outgoing halfedges. This _should_ touch every halfedge.
      size_t currHe = vHalfedge[iV];
      size_t firstHe = currHe;
      do {

        GC_SAFETY_ASSERT(!halfedgeSeen[currHe], "somehow encountered outgoing halfedge before orbiting v");
        halfedgeSeen[currHe] = true;

        currHe = heNext[heTwin(currHe)];
      } while (currHe != firstHe);
    }

    // Verify that we actually did touch every halfedge.
    for (size_t iHe = 0; iHe < nHalfedgesCount; iHe++) {
      GC_SAFETY_ASSERT(halfedgeSeen[currHe], "mesh not manifold. Vertex " + std::to_string(heVertex[iHe]) +
                                                 " has disconnected neighborhoods incident (imagine an hourglass)");
    }
  }
#endif


  // Print some nice statistics
  if (verbose) {
    std::cout << "Constructed halfedge mesh with: " << std::endl;
    std::cout << "    # verts =  " << nVertices() << std::endl;
    std::cout << "    # edges =  " << nEdges() << std::endl;
    std::cout << "    # faces =  " << nFaces() << std::endl;
    std::cout << "    # halfedges =  " << nHalfedges() << std::endl;
    std::cout << "      and " << nBoundaryLoops() << " boundary components. " << std::endl;
    std::cout << "Construction took " << pretty_time(FINISH_TIMING(construction)) << std::endl;
  }
}


HalfedgeMesh::HalfedgeMesh(const std::vector<std::vector<size_t>>& polygons) {
  /*   High-level outline of this algorithm:
   *
   *      0. Count how many of each object we will need so we can pre-allocate them
   *
   *      1. Iterate over the faces of the input mesh, creating edge, face, and halfedge objects
   *
   *      2. Walk around boundaries, marking boundary edge and creating imaginary halfedges/faces
   *
   *      3. Copy the vertex positions to the geometry associated with the new mesh
   *
   */

  START_TIMING(construction)

  // === 0. Count needed objects

  // Find which vertices are actually used
  size_t maxVertexID = 0;
  for (auto poly : polygons) {
    for (auto i : poly) {
      maxVertexID = std::max(maxVertexID, i);
    }
  }
  std::vector<bool> usedVerts(maxVertexID);
  std::vector<size_t> usedVertsIndex(maxVertexID);
  std::fill(usedVerts.begin(), usedVerts.end(), false);
  std::fill(usedVertsIndex.begin(), usedVertsIndex.end(), 0);
  size_t nVerts = 0;
  for (auto poly : polygons) {
    for (auto i : poly) {
      usedVerts[i] = true;
    }
  }
  for (size_t i = 0; i < maxVertexID; i++) {
    if (usedVerts[i]) {
      usedVertsIndex[i] = nVerts++;
    }
  }

  // === 0.5 Efficiently build an adjacent-face lookup table

  size_t INVALID_IND = std::numeric_limits<size_t>::max();

  // Build a sorted list of (directed halfedge) neighbors of each vertex in compressed format.

  // Count neighbors of each vertex
  size_t nDirected = 0;
  size_t maxPolyDegree = 0;
  std::vector<size_t> vertexNeighborsCount(maxVertexID, 0);
  std::vector<size_t> vertexNeighborsStart(maxVertexID + 1);
  for (size_t iFace = 0; iFace < polygons.size(); iFace++) {
    auto poly = polygons[iFace];
    nDirected += poly.size();
    maxPolyDegree = std::max(maxPolyDegree, poly.size());
    for (size_t j : poly) {
      vertexNeighborsCount[j]++;
    }
  }

  // Build a running sum of the number of neighbors to use a compressed list
  vertexNeighborsStart[0] = 0;
  size_t runningSum = 0;
  for (size_t iVert = 0; iVert < maxVertexID; iVert++) {
    runningSum += vertexNeighborsCount[iVert];
    vertexNeighborsStart[iVert + 1] = runningSum;
  }

  // Populate the compressed list
  // Each vertex's neighbors are stored between vertexNeighborsStart[i] and
  // vertexNeighborsStart[i+1]
  std::vector<size_t> vertexNeighborsInd(vertexNeighborsStart.begin(), vertexNeighborsStart.end() - 1);
  std::vector<size_t> allVertexNeighbors(nDirected);
  for (size_t iFace = 0; iFace < polygons.size(); iFace++) {
    auto poly = polygons[iFace];
    for (size_t j = 0; j < poly.size(); j++) {
      size_t fromInd = poly[j];
      size_t toInd = poly[(j + 1) % poly.size()];
      allVertexNeighbors[vertexNeighborsInd[fromInd]] = toInd;
      vertexNeighborsInd[fromInd]++;
    }
  }

  // Sort each of the sublists in the compressed list
  for (size_t iVert = 0; iVert < maxVertexID; iVert++) {
    std::sort(allVertexNeighbors.begin() + vertexNeighborsStart[iVert],
              allVertexNeighbors.begin() + vertexNeighborsStart[iVert + 1]);
  }

  // Count real and imaginary edges and faces and cache adjacent twin indices
  // Note: counting boundary loops is kinda difficult, so we wait to do so until
  // the final step
  size_t nPairedEdges = 0;
  size_t nUnpairedEdges = 0;
  std::vector<size_t> twinInd(allVertexNeighbors.size());
  for (size_t iVert = 0; iVert < maxVertexID; iVert++) {
    size_t jStart = vertexNeighborsStart[iVert];
    size_t jEnd = vertexNeighborsStart[iVert + 1];
    for (size_t jInd = jStart; jInd < jEnd; jInd++) {
      size_t jVert = allVertexNeighbors[jInd];

      // Search for the j --> i edge
      size_t searchResult =
          halfedgeLookup(allVertexNeighbors, iVert, vertexNeighborsStart[jVert], vertexNeighborsStart[jVert + 1]);
      twinInd[jInd] = searchResult;
      if (searchResult == INVALID_IND) {
        nUnpairedEdges++;
      } else {
        nPairedEdges++;
      }
    }
  }
  nPairedEdges /= 2;

  size_t nTotalEdges = nPairedEdges + nUnpairedEdges;
  size_t nRealHalfedgesToMake =
      2 * nPairedEdges + nUnpairedEdges; // use a temp variable here because this is tracked as part of class state
  size_t nImaginaryHalfedges = nUnpairedEdges;
  size_t nRealFaces = polygons.size();

  // Allocate space and construct elements
  rawHalfedges.reserve(nRealHalfedgesToMake + nImaginaryHalfedges);
  for (size_t i = 0; i < nRealHalfedgesToMake; i++) getNewHalfedge(true);
  for (size_t i = 0; i < nImaginaryHalfedges; i++) getNewHalfedge(false);
  rawVertices.reserve(nVerts);
  for (size_t i = 0; i < nVerts; i++) getNewVertex();
  rawEdges.reserve(nTotalEdges);
  for (size_t i = 0; i < nTotalEdges; i++) getNewEdge();
  rawFaces.reserve(nRealFaces);
  for (size_t i = 0; i < nRealFaces; i++) getNewFace();

  // === 1. Create faces, edges, and halfedges

  // Keep track of the edges we've already created since we only need one per
  // edge
  std::vector<Edge> sharedEdges(twinInd.size(), Edge());

  // Iterate over faces
  size_t iFace = 0;
  size_t iEdge = 0;
  size_t iHalfedge = 0;
  std::vector<Halfedge> thisFaceHalfedges;
  thisFaceHalfedges.reserve(maxPolyDegree);
  for (auto poly : polygons) {
    size_t degree = poly.size();

    // Create a new face object
    Face f{&rawFaces[iFace]};
    f->isReal = true;

    // The halfedges that make up this face
    // std::vector<Halfedge> thisFaceHalfedges(degree);
    thisFaceHalfedges.resize(degree);

    for (size_t iPolyEdge = 0; iPolyEdge < degree; iPolyEdge++) {
      Halfedge he{&rawHalfedges[iHalfedge]};
      iHalfedge++;

      size_t ind1 = poly[iPolyEdge];
      size_t ind2 = poly[(iPolyEdge + 1) % degree];

      // Connect up pointers
      he->vertex = &rawVertices[usedVertsIndex[ind1]];
      he->vertex->halfedge = he.ptr;
      he->face = f.ptr;
      thisFaceHalfedges[iPolyEdge] = he;
      he->twin = nullptr; // ensure this is null so we can detect boundaries below

      // Get a reference to the edge shared by this and its twin, creating the
      // object if needed
      size_t myHeInd =
          halfedgeLookup(allVertexNeighbors, ind2, vertexNeighborsStart[ind1], vertexNeighborsStart[ind1 + 1]);
      size_t twinHeInd = twinInd[myHeInd];
      bool edgeAlreadyCreated = (twinHeInd != INVALID_IND) && (sharedEdges[twinHeInd] != Edge());

      if (edgeAlreadyCreated) {
        Edge sharedEdge = sharedEdges[twinHeInd];
        he->edge = sharedEdge.ptr;
        he->twin = sharedEdge->halfedge;
        he->twin->twin = he.ptr;
      } else {
        Edge sharedEdge{&rawEdges[iEdge]};
        iEdge++;
        sharedEdge->halfedge = he.ptr;
        he->edge = sharedEdge.ptr;
        sharedEdges[myHeInd] = sharedEdge;
      }
    }

    // Do one more lap around the face to set next pointers
    for (size_t iPolyEdge = 0; iPolyEdge < degree; iPolyEdge++) {
      thisFaceHalfedges[iPolyEdge]->next = thisFaceHalfedges[(iPolyEdge + 1) % degree].ptr;
    }

    f->halfedge = thisFaceHalfedges[0].ptr;
    thisFaceHalfedges.clear();
    iFace++;
  }

  // === 2. Walk the boundary to find/create boundary cycles

  // First, do a pre-walk to count the boundary loops we will need and allocate
  // them
  size_t nBoundaryLoops = 0;
  std::set<Halfedge> walkedHalfedges;
  for (size_t iHe = 0; iHe < nRealHalfedges(); iHe++) {
    if (halfedge(iHe)->twin == nullptr && walkedHalfedges.find(halfedge(iHe)) == walkedHalfedges.end()) {
      nBoundaryLoops++;
      Halfedge currHe = halfedge(iHe);
      walkedHalfedges.insert(currHe);
      size_t walkCount = 0;
      do {
        currHe = currHe->next;
        while (currHe->twin != nullptr) {
          currHe = currHe->twin->next;
          walkCount++;
          if (walkCount > nRealHalfedges()) {
            throw std::runtime_error(
                "Encountered infinite loop while constructing halfedge mesh. Are you sure the input is manifold?");
          }
        }
        walkedHalfedges.insert(currHe);
      } while (currHe != halfedge(iHe));
    }
  }
  rawBoundaryLoops.resize(nBoundaryLoops);


  // Now do the actual walk in which we construct and connect objects
  size_t iBoundaryLoop = 0;
  size_t iImaginaryHalfedge = 0;
  std::vector<char> vertexAppearedInBoundaryLoop(nVerts, false);
  for (size_t iHe = 0; iHe < nRealHalfedges(); iHe++) {
    // Note: If distinct holes share a given vertex, this algorithm will see
    // them as a single "figure 8-like"
    // hole and connect them with a single imaginary face.
    // TODO: fix this?

    // If this halfedge doesn't have a twin, it must be on a boundary (or have
    // already been processed while walking a hole)
    if (halfedge(iHe)->twin == nullptr) {
      // Create a boundary loop for this hole
      BoundaryLoop boundaryLoop{&rawBoundaryLoops[iBoundaryLoop]};
      boundaryLoop->isReal = false;

      // Walk around the boundary loop, creating imaginary halfedges
      Halfedge currHe = halfedge(iHe);
      Halfedge prevHe{nullptr};
      bool finished = false;
      while (!finished) {

        // Check for non-manifoldness via non-disk-like boundary vertex
        size_t iV = currHe.vertex().ptr - &rawVertices[0];
        if (vertexAppearedInBoundaryLoop[iV] == true) {
          throw std::runtime_error("Input mesh is nonmanifold: vertex appears in two distinct boundary loops");
        }
        vertexAppearedInBoundaryLoop[iV] = true;

        // Create a new, imaginary halfedge
        Halfedge newHe{&rawHalfedges[nRealHalfedgesCount + iImaginaryHalfedge]};
        boundaryLoop->halfedge = newHe.ptr;
        iImaginaryHalfedge++;

        // Connect up pointers
        newHe->isReal = false;
        newHe->twin = currHe.ptr;
        currHe->twin = newHe.ptr;
        newHe->face = boundaryLoop.ptr;
        newHe->edge = currHe->edge;
        newHe->vertex = currHe->next->vertex;
        currHe->vertex->halfedge =
            currHe.ptr; // ensure that halfedge for boundary vertex is the one that starts the boundary

        // Some pointers need values only visible from the previous iteration of
        // the loop.
        // The first one we process gets missed, handle it at the end outside
        // the loop.
        if (prevHe == nullptr) {
          newHe->next = nullptr;
        } else {
          newHe->next = prevHe->twin;
        }

        // Set the isBoundary property where appropriate
        currHe->face->isBoundary = true;
        currHe->vertex->isBoundary = true;
        currHe->edge->isBoundary = true;

        // Prepare for the next iteration
        prevHe = currHe;
        currHe = currHe->next;
        while (currHe->twin != nullptr) {
          // When we've finished walking the loop, we'll be able to tell
          // because we spin in circles around the vertex trying to continue
          // the loop. Detect that here and quit.
          if (currHe->twin->next == nullptr) {
            finished = true;
            break;
          }

          currHe = currHe->twin->next;
        }
      }

      // As noted above, the pointers don't get set properly on the first
      // iteration of
      // the loop above because we don't have a reference to prev yet. Fix that
      // here.
      halfedge(iHe)->twin->next = prevHe->twin;

      iBoundaryLoop++;
    }
  }

  { // Check that the input was manifold in the sense that each vertex has a single connected loop of faces around it.
    // This is just a sanity check, and can be skipped if the input is trusted. However, if not checked, we could output
    // a nonmanifold "mesh".
    std::vector<char> vertexSeen(nVerts, false);
    std::vector<char> halfedgeSeen(nRealHalfedges() + this->nImaginaryHalfedges(), false);
    for (size_t iHe = 0; iHe < nRealHalfedges(); iHe++) {
      if (halfedgeSeen[iHe]) continue;

      Halfedge* firstHe = &rawHalfedges[iHe];
      Halfedge* currHe = firstHe;
      size_t iV = firstHe->vertex - &rawVertices[0];
      if (vertexSeen[iV]) {
        throw std::runtime_error(
            "Vertex neighborhood is nonmanifold: >1 distinct ring of triangles incident on vertex.");
      }
      vertexSeen[iV] = true;
      do {
        size_t currIHe = currHe - &rawHalfedges[0];
        halfedgeSeen[currIHe] = true;
        currHe = currHe->twin->next;
      } while (currHe != firstHe);
    }
  }


// When in debug mode, mesh elements know what mesh they are a part of so
// we can do assertions for saftey checks.
#ifndef NDEBUG
  for (Halfedge x : allHalfedges()) {
    x->parentMesh = this;
  }
  for (Vertex x : vertices()) {
    x->parentMesh = this;
  }
  for (Edge x : edges()) {
    x->parentMesh = this;
  }
  for (Face x : faces()) {
    x->parentMesh = this;
  }
  for (Face x : boundaryLoops()) {
    x->parentMesh = this;
  }
#endif

  // === 3. Map vertices in the halfedge mesh to the associated vertex
  // coordinates in space
  // Create the vertex objects and build a map to find them
  size_t iVert = 0;
  geometry = new Geometry<Euclidean>(*this);
  for (size_t i = 0; i < maxVertexID; i++) {
    if (usedVerts[i]) {
      geometry->position(vertex(iVert)) = input.vertexCoordinates[i];
      iVert++;
    }
  }

  // Print some nice statistics
  std::cout << "Constructed halfedge mesh with: " << std::endl;
  std::cout << "    # verts =  " << nVertices() << std::endl;
  std::cout << "    # edges =  " << nEdges() << std::endl;
  std::cout << "    # faces =  " << nFaces() << std::endl;
  std::cout << "    # halfedges =  " << nRealHalfedges() + this->nImaginaryHalfedges() << std::endl;
  std::cout << "      and " << nBoundaryLoops << " boundary components. " << std::endl;
  std::cout << "Construction took " << pretty_time(FINISH_TIMING(construction)) << std::endl;


  // Compute some basic information about the mesh
}

HalfedgeMesh::~HalfedgeMesh() {
  for (auto& f : meshDeleteCallbackList) {
    f();
  }
}


// Returns true if and only if all faces are triangles
bool HalfedgeMesh::isSimplicial() {
  for (Face f : faces()) {
    if (f.degree() != 3) {
      return false;
    }
  }
  return true;
}

size_t HalfedgeMesh::nConnectedComponents() {
  VertexData<size_t> vertInd = getVertexIndices();
  DisjointSets dj(nVertices());
  for (Edge e : edges()) {
    dj.merge(vertInd[e.halfedge().vertex()], vertInd[e.halfedge().twin().vertex()]);
  }
  std::unordered_set<size_t> distinctComponents;
  for (size_t i = 0; i < nVertices(); i++) {
    distinctComponents.insert(dj.find(i));
  }
  return distinctComponents.size();
}

size_t HalfedgeMesh::nInteriorVertices() {
  size_t _nInteriorVertices = 0;
  for (const Vertex v : vertices()) {
    if (!v->isBoundary) {
      _nInteriorVertices++;
    }
  }
  return _nInteriorVertices;
}

// Counts the total number of faces in a triangulation of this mesh,
// corresponding to the triangulation given by Face::triangulate().
size_t HalfedgeMesh::nFacesTriangulation() {
  size_t _nFacesTriangulation = 0;
  for (Face f : faces()) {
    _nFacesTriangulation += f.degree() - 2;
  }
  return _nFacesTriangulation;
}


int HalfedgeMesh::eulerCharacteristic() {
  // be sure to do intermediate arithmetic with large, signed integers
  return static_cast<int>(static_cast<long long int>(nVertices()) - static_cast<long long int>(nEdges()) +
                          static_cast<long long int>(nFaces()));
}

size_t HalfedgeMesh::longestBoundaryLoop() {
  int max = 0;
  size_t _longestBoundaryLoop = 0;
  for (size_t i = 0; i < nBoundaryLoops(); i++) {
    // Count halfedges on boundary and check if greater than max
    int n = 0;
    Halfedge he = boundaryLoop(i).halfedge();
    Halfedge h = he;
    do {
      n++;
      h = h.next();
    } while (h != he);

    if (n > max) {
      max = n;
      _longestBoundaryLoop = i;
    }
  }
  return _longestBoundaryLoop;
}

VertexData<size_t> HalfedgeMesh::getVertexIndices() {
  VertexData<size_t> indices(this);
  size_t i = 0;
  for (Vertex v : vertices()) {
    indices[v] = i;
    i++;
  }
  return indices;
}

VertexData<size_t> HalfedgeMesh::getInteriorVertexIndices() {
  VertexData<size_t> indices(this);
  size_t i = 0;
  for (Vertex v : vertices()) {
    if (v->isBoundary) {
      indices[v] = -7;
    } else {
      indices[v] = i;
      i++;
    }
  }
  return indices;
}

FaceData<size_t> HalfedgeMesh::getFaceIndices() {
  FaceData<size_t> indices(this);
  size_t i = 0;
  for (Face f : faces()) {
    indices[f] = i;
    i++;
  }
  return indices;
}

EdgeData<size_t> HalfedgeMesh::getEdgeIndices() {
  EdgeData<size_t> indices(this);
  size_t i = 0;
  for (Edge e : edges()) {
    indices[e] = i;
    i++;
  }
  return indices;
}

HalfedgeData<size_t> HalfedgeMesh::getHalfedgeIndices() {
  HalfedgeData<size_t> indices(this);
  size_t i = 0;
  for (Halfedge he : allHalfedges()) {
    indices[he] = i;
    i++;
  }
  return indices;
}

CornerData<size_t> HalfedgeMesh::getCornerIndices() {
  CornerData<size_t> indices(this);
  size_t i = 0;
  for (Corner c : corners()) {
    if (c.halfedge().isReal()) {
      indices[c] = i;
      i++;
    }
  }
  return indices;
}

HalfedgeMesh* HalfedgeMesh::copy() {

  HalfedgeMeshDataTransfer t;
  return copy(t);
}


HalfedgeMesh* HalfedgeMesh::copy(HalfedgeMeshDataTransfer& dataTransfer) {

  HalfedgeMesh* newMesh = new HalfedgeMesh();

  // Copy vectors
  newMesh->rawHalfedges = rawHalfedges;
  newMesh->rawVertices = rawVertices;
  newMesh->rawEdges = rawEdges;
  newMesh->rawFaces = rawFaces;
  newMesh->rawBoundaryLoops = rawBoundaryLoops;

  newMesh->nRealHalfedgesCount = nRealHalfedgesCount;
  newMesh->nImaginaryHalfedgesCount = nImaginaryHalfedgesCount;
  newMesh->nVerticesCount = nVerticesCount;
  newMesh->nEdgesCount = nEdgesCount;
  newMesh->nFacesCount = nFacesCount;
  newMesh->nBoundaryLoopsCount = nBoundaryLoopsCount;

  // Copy other values
  newMesh->nextElemID = nextElemID;
  newMesh->isCompressedFlag = isCompressedFlag;
  newMesh->isCanonicalFlag = isCanonicalFlag;

  // Create the data transfer object
  dataTransfer = HalfedgeMeshDataTransfer(this, newMesh);

  // Build maps
  for (size_t i = 0; i < rawHalfedges.size(); i++) {
    dataTransfer.heMap[halfedge(i)] = &newMesh->rawHalfedges[i];
  }
  for (size_t i = 0; i < rawVertices.size(); i++) {
    dataTransfer.vMap[vertex(i)] = &newMesh->rawVertices[i];
  }
  for (size_t i = 0; i < rawEdges.size(); i++) {
    dataTransfer.eMap[edge(i)] = &newMesh->rawEdges[i];
  }
  for (size_t i = 0; i < rawFaces.size(); i++) {
    dataTransfer.fMap[face(i)] = &newMesh->rawFaces[i];
  }
  for (size_t i = 0; i < rawBoundaryLoops.size(); i++) {
    dataTransfer.fMap[boundaryLoop(i)] = &newMesh->rawBoundaryLoops[i];
  }

  // Shift pointers
  for (size_t i = 0; i < rawHalfedges.size(); i++) {
    newMesh->rawHalfedges[i].next = dataTransfer.heMap[halfedge(i).next()].ptr;
    newMesh->rawHalfedges[i].twin = dataTransfer.heMap[halfedge(i).twin()].ptr;
    newMesh->rawHalfedges[i].vertex = dataTransfer.vMap[halfedge(i).vertex()].ptr;
    newMesh->rawHalfedges[i].edge = dataTransfer.eMap[halfedge(i).edge()].ptr;
    newMesh->rawHalfedges[i].face = dataTransfer.fMap[halfedge(i).face()].ptr;
  }
  for (size_t i = 0; i < rawVertices.size(); i++) {
    newMesh->rawVertices[i].halfedge = dataTransfer.heMap[vertex(i).halfedge()].ptr;
  }
  for (size_t i = 0; i < rawEdges.size(); i++) {
    newMesh->rawEdges[i].halfedge = dataTransfer.heMap[edge(i).halfedge()].ptr;
  }
  for (size_t i = 0; i < rawFaces.size(); i++) {
    newMesh->rawFaces[i].halfedge = dataTransfer.heMap[face(i).halfedge()].ptr;
  }
  for (size_t i = 0; i < rawBoundaryLoops.size(); i++) {
    newMesh->rawBoundaryLoops[i].halfedge = dataTransfer.heMap[boundaryLoop(i).halfedge()].ptr;
  }


#ifndef NDEBUG
  // Set the parent mesh pointers to point to the correct mesh
  for (size_t i = 0; i < newMesh->rawHalfedges.size(); i++) {
    newMesh->rawHalfedges[i].parentMesh = newMesh;
  }
  for (size_t i = 0; i < newMesh->rawVertices.size(); i++) {
    newMesh->rawVertices[i].parentMesh = newMesh;
  }
  for (size_t i = 0; i < newMesh->rawEdges.size(); i++) {
    newMesh->rawEdges[i].parentMesh = newMesh;
  }
  for (size_t i = 0; i < newMesh->rawFaces.size(); i++) {
    newMesh->rawFaces[i].parentMesh = newMesh;
  }
  for (size_t i = 0; i < newMesh->rawBoundaryLoops.size(); i++) {
    newMesh->rawBoundaryLoops[i].parentMesh = newMesh;
  }
#endif

  dataTransfer.generateReverseMaps();

  return newMesh;
}

std::vector<std::vector<size_t>> HalfedgeMesh::getPolygonSoupFaces() {

  std::vector<std::vector<size_t>> result;

  VertexData<size_t> vInd = getVertexIndices();
  for (Face f : faces()) {
    std::vector<size_t> faceList;
    for (Vertex v : f.adjacentVertices()) {
      faceList.push_back(vInd[v]);
    }
    result.push_back(faceList);
  }

  return result;
}

bool HalfedgeMesh::flip(Edge eFlip) {

  Edge* e = eFlip.ptr;

  // Get halfedges of first face
  Halfedge* ha1 = e->halfedge;
  if (!ha1->isReal) return false; // don't flip boundary edges
  Halfedge* ha2 = ha1->next;
  Halfedge* ha3 = ha2->next;
  if (ha3->next != ha1) return false; // not a triangle

  // Get halfedges of second face
  Halfedge* hb1 = ha1->twin;
  if (!hb1->isReal) return false; // don't flip boundary edges
  Halfedge* hb2 = hb1->next;
  Halfedge* hb3 = hb2->next;
  if (hb3->next != hb1) return false; // not a triangle

  if (ha2 == hb1 || hb2 == ha1) return false; // incident on degree 1 vertex

  // Get vertices and faces
  Vertex* va = ha1->vertex;
  Vertex* vb = hb1->vertex;
  Vertex* vc = ha3->vertex;
  Vertex* vd = hb3->vertex;
  Face* fa = ha1->face;
  Face* fb = hb1->face;

  // Update vertex pointers
  if (va->halfedge == ha1) va->halfedge = hb2;
  if (vb->halfedge == hb1) vb->halfedge = ha2;
  // (vc and vd can't be invalidated by the flip)

  // Update edge pointers
  // (e still has the same halfedges)

  // Update face pointers
  fa->halfedge = ha1;
  fb->halfedge = hb1;

  // Update halfedge pointers
  ha1->next = hb3;
  hb3->next = ha2;
  ha2->next = ha1;
  hb1->next = ha3;
  ha3->next = hb2;
  hb2->next = hb1;
  ha1->vertex = vc;
  hb1->vertex = vd;
  ha3->face = fb;
  hb3->face = fa;

  isCanonicalFlag = false;
  return true;
}

Halfedge HalfedgeMesh::insertVertexAlongEdge(Edge eIn) {

  DynamicEdge eInD(eIn, this);

  // == Gather / create elements
  // Faces are identified as 'A', and 'B'

  // Create first, because getNew() could invalidate pointers
  Vertex* newV = getNewVertex();
  Edge* newE = getNewEdge();
  Halfedge* heANew = getNewHalfedge(true);
  Halfedge* heBNew = getNewHalfedge(true);

  Edge e = eInD;
  Halfedge* heACenter = e.halfedge().ptr;
  Halfedge* heBCenter = heACenter->twin;
  // Halfedge* heANext = heACenter->next;
  Halfedge* heBNext = heBCenter->next;
  Halfedge* heAPrev = Halfedge{heACenter}.prev().ptr;
  // Halfedge* heBPrev = Halfedge{heBCenter}.prev().ptr;
  Face* fA = heACenter->face;
  Face* fB = heBCenter->face;
  Vertex* oldVBottom = heACenter->vertex;

  // == Hook up all the pointers

  // New vertex
  newV->halfedge = heACenter;

  // New edge
  newE->halfedge = heANew;

  // New halfedge A
  heANew->twin = heBNew;
  heANew->next = heACenter;
  heANew->vertex = oldVBottom;
  heANew->edge = newE;
  heANew->face = fA;

  // New halfedge B
  heBNew->twin = heANew;
  heBNew->next = heBNext;
  heBNew->vertex = newV;
  heBNew->edge = newE;
  heBNew->face = fB;

  // Fix pointers for old halfedges
  heBCenter->next = heBNew;
  heAPrev->next = heANew;
  oldVBottom->halfedge = heBNext;
  heACenter->vertex = newV;

  isCanonicalFlag = false;

  return Halfedge{heACenter};
}


Vertex HalfedgeMesh::splitEdge(Edge e) {

  // Validate that faces are triangular and real
  if (e.isBoundary() || e.halfedge().face().degree() != 3 || e.halfedge().twin().face().degree() != 3 ||
      e.halfedge().face() == e.halfedge().twin().face()) {
    throw std::logic_error("Can only split non-boundary edge which borders two distinct triangular faces");
  }

  // First operation: insert a new vertex along the edge
  Vertex newV = insertVertexAlongEdge(e).vertex();


  // Second operation: connect both of the new faces
  DynamicFace fOppA(newV.halfedge().face(), this);
  DynamicVertex vOppA(newV.halfedge().next().next().vertex(), this);
  DynamicFace fOppB(newV.halfedge().twin().face(), this);
  DynamicVertex vOppB(newV.halfedge().twin().next().next().next().vertex(), this);

  connectVertices(fOppA, vOppA, newV);
  connectVertices(fOppB, vOppB, newV);

  isCanonicalFlag = false;
  return newV;
}

Halfedge HalfedgeMesh::splitEdgeReturnHalfedge(Edge e) {

  // Validate that faces are triangular and real
  if (e.isBoundary() || e.halfedge().face().degree() != 3 || e.halfedge().twin().face().degree() != 3 ||
      e.halfedge().face() == e.halfedge().twin().face()) {
    throw std::logic_error("Can only split non-boundary edge which borders two distinct triangular faces");
    // note: if removing boundary restriction, need to fix return value below
  }

  // Save this
  DynamicHalfedge heIn(e.halfedge(), this);

  // First operation: insert a new vertex along the edge
  Vertex newV = insertVertexAlongEdge(e).vertex();


  // Second operation: connect both of the new faces
  DynamicFace fOppA(newV.halfedge().face(), this);
  DynamicVertex vOppA(newV.halfedge().next().next().vertex(), this);
  DynamicFace fOppB(newV.halfedge().twin().face(), this);
  DynamicVertex vOppB(newV.halfedge().twin().next().next().next().vertex(), this);

  connectVertices(fOppA, vOppA, newV);
  connectVertices(fOppB, vOppB, newV);

  // Find the useful halfedge pointer to return using from the halfedge we saved
  Halfedge toReturn = Halfedge(heIn).twin().next().twin().next().twin();

  isCanonicalFlag = false;
  return toReturn;
}


Halfedge HalfedgeMesh::connectVertices(Vertex vA, Vertex vB) {

  // Find the shared face and call the main version
  std::unordered_set<Face> aFaces;
  for (Face f : vA.adjacentFaces()) {
    aFaces.insert(f);
  }
  Face sharedFace = Face();
  for (Face f : vB.adjacentFaces()) {
    if (aFaces.find(f) != aFaces.end()) {
      sharedFace = f;
      break;
    }
  }

  if (sharedFace == Face()) {
    throw std::logic_error("Vertices do not contain shared face");
  }

  isCanonicalFlag = false;
  return connectVertices(sharedFace, vA, vB);
}

Halfedge HalfedgeMesh::tryConnectVertices(Vertex vA, Vertex vB) {


  // TODO much of the connectVertices() logic is O(N_VERTICES_IN_FACE) even though it doesn't really need to be.

  // Early-out if same
  if (vA == vB) {
    return Halfedge();
  }

  // Find the shared face and call the main version
  std::unordered_set<Face> aFaces;
  for (Face f : vA.adjacentFaces()) {
    aFaces.insert(f);
  }
  Face sharedFace = Face();
  for (Face f : vB.adjacentFaces()) {
    if (aFaces.find(f) != aFaces.end()) {
      sharedFace = f;
      break;
    }
  }

  // Fail if no shared face
  if (sharedFace == Face()) {
    return Halfedge();
  }

  // Check if adjacent
  for (Halfedge he : sharedFace.adjacentHalfedges()) {
    if ((he.vertex() == vA && he.twin().vertex() == vB) || (he.vertex() == vB && he.twin().vertex() == vA)) {
      return Halfedge();
    }
  }

  isCanonicalFlag = false;
  return connectVertices(sharedFace, vA, vB);
}

Halfedge HalfedgeMesh::tryConnectVertices(Vertex vA, Vertex vB, Face face) {


  // TODO much of the connectVertices() logic is O(N_VERTICES_IN_FACE) even though it doesn't really need to be.

  // Early-out if same
  if (vA == vB) {
    return Halfedge();
  }

  // Find the shared face and call the main version

  // Check if adjacent
  bool foundA = false;
  bool foundB = false;
  for (Halfedge he : face.adjacentHalfedges()) {
    if ((he.vertex() == vA && he.twin().vertex() == vB) || (he.vertex() == vB && he.twin().vertex() == vA)) {
      return Halfedge();
    }

    if (he.vertex() == vA) {
      foundA = true;
    }
    if (he.vertex() == vB) {
      foundB = true;
    }
  }

  // One of the vertices isn't in the face we're supposed to work in
  if (!foundA) {
    return Halfedge();
  }
  if (!foundB) {
    return Halfedge();
  }

  isCanonicalFlag = false;
  return connectVertices(face, vA, vB);
}

Halfedge HalfedgeMesh::connectVertices(Face faceIn, Vertex vAIn, Vertex vBIn) {

  DynamicFace faceInD(faceIn, this);
  DynamicVertex vAInD(vAIn, this);
  DynamicVertex vBInD(vBIn, this);

  // == Create new elements
  Halfedge* heANew = getNewHalfedge(true);
  Halfedge* heBNew = getNewHalfedge(true);
  Edge* eNew = getNewEdge();
  Face* fB = getNewFace();

  Face face(faceInD);
  Vertex vA(vAInD);
  Vertex vB(vBInD);

  // == Find useful halfedges around the face
  Halfedge* heANext;
  Halfedge* heBNext;
  Halfedge* heAPrev;
  Halfedge* heBPrev;
  for (Halfedge he : face.adjacentHalfedges()) {
    if (he.vertex() == vA) {
      heANext = he.ptr;
    }
    if (he.vertex() == vB) {
      heBNext = he.ptr;
    }
  }
  for (Halfedge he : face.adjacentHalfedges()) {
    if (he.next().ptr == heANext) {
      heAPrev = he.ptr;
    }
    if (he.next().ptr == heBNext) {
      heBPrev = he.ptr;
    }
  }

  // == Detect bad cases
  if (vA == vB) throw std::logic_error("Tried to connect vertex to self");
  if (heANext == heBPrev || heBNext == heBPrev) throw std::logic_error("Tried to connect adjacent vertices");

  // == Gather other elements
  Face* fA = heBNext->face;
  Vertex* vAp = vA.ptr;
  Vertex* vBp = vB.ptr;

  // == Hook up all the pointers

  // Faces
  fA->halfedge = heBNext;
  fB->halfedge = heANext;

  // Vertices
  vAp->halfedge = heANew;
  vBp->halfedge = heBNew;

  // New edge
  eNew->halfedge = heANew;

  // Halfedges
  heANew->twin = heBNew;
  heANew->next = heBNext;
  heANew->vertex = vAp;
  heANew->edge = eNew;
  heANew->face = fA;

  heBNew->twin = heANew;
  heBNew->next = heANext;
  heBNew->vertex = vBp;
  heBNew->edge = eNew;
  heBNew->face = fB;

  heAPrev->next = heANew;
  heBPrev->next = heBNew;

  // Set all other new .face pointers to fB
  Halfedge* currHe = heANext;
  while (currHe != heBNew) {
    currHe->face = fB;
    currHe = currHe->next;
  }

  isCanonicalFlag = false;
  return heANew;
}


Vertex HalfedgeMesh::insertVertex(Face fIn) {

  DynamicFace fInD(fIn, this);

  // Create the new center vertex
  DynamicVertex centerVertD(getNewVertex(), this);


  // Count degree to allocate elements
  size_t faceDegree = 0;
  for (Halfedge he : Face(fInD).adjacentHalfedges()) {
    faceDegree++;
  }

  // == Create new halfedges/edges/faces around the center vertex

  // Create all of the new elements first, then hook them up below, as this can invalidate pointers.
  std::vector<DynamicFace> innerFacesD;
  std::vector<DynamicHalfedge> leadingHalfedgesD; // the one that points towards the center
  std::vector<DynamicHalfedge> trailingHalfedgesD;
  std::vector<DynamicEdge> innerEdgesD; // aligned with leading he
  for (size_t i = 0; i < faceDegree; i++) {
    // Re-use first face
    if (i == 0) {
      innerFacesD.push_back(fInD);
    } else {
      innerFacesD.push_back(DynamicFace(getNewFace(), this));
    }

    leadingHalfedgesD.push_back(DynamicHalfedge(getNewHalfedge(true), this));
    trailingHalfedgesD.push_back(DynamicHalfedge(getNewHalfedge(true), this));
    innerEdgesD.push_back(DynamicEdge(getNewEdge(), this));
  }

  // Copy to raw pointers now that they are allocated
  std::vector<Face*> innerFaces;
  std::vector<Halfedge*> leadingHalfedges; // the one that points towards the center
  std::vector<Halfedge*> trailingHalfedges;
  std::vector<Edge*> innerEdges; // aligned with leading he
  for (size_t i = 0; i < faceDegree; i++) {
    innerFaces.push_back(Face(innerFacesD[i]).ptr);
    leadingHalfedges.push_back(Halfedge(leadingHalfedgesD[i]).ptr);
    trailingHalfedges.push_back(Halfedge(trailingHalfedgesD[i]).ptr);
    innerEdges.push_back(Edge(innerEdgesD[i]).ptr);
  }
  Vertex centerVert = centerVertD;

  // Form this list before we start, because we're about to start breaking pointers
  std::vector<Halfedge*> faceBoundaryHalfedges;
  for (Halfedge he : Face(fInD).adjacentHalfedges()) {
    faceBoundaryHalfedges.push_back(he.ptr);
  }

  // Connect up all the pointers
  // Each iteration processes one inner face
  for (size_t i = 0; i < faceDegree; i++) {

    // Gather pointers
    Face* f = innerFaces[i];
    Edge* e = innerEdges[i];
    Edge* prevE = innerEdges[(i + faceDegree - 1) % faceDegree];
    Halfedge* leadingHe = leadingHalfedges[i];
    Halfedge* trailingHe = trailingHalfedges[i];
    Halfedge* boundaryHe = faceBoundaryHalfedges[i];
    Halfedge* nextTrailingHe = trailingHalfedges[(i + 1) % faceDegree];
    Halfedge* prevLeadingHe = leadingHalfedges[(i + faceDegree - 1) % faceDegree];

    // face
    f->halfedge = boundaryHe;
    f->isReal = true;

    // edge
    e->halfedge = leadingHe;

    // leading halfedge
    leadingHe->twin = nextTrailingHe;
    leadingHe->next = trailingHe;
    leadingHe->vertex = boundaryHe->next->vertex;
    leadingHe->edge = e;
    leadingHe->face = f;

    // trailing halfedge
    trailingHe->twin = prevLeadingHe;
    trailingHe->next = boundaryHe;
    trailingHe->vertex = centerVert.ptr;
    trailingHe->edge = prevE;
    trailingHe->face = f;

    // boundary halfedge
    boundaryHe->next = leadingHe;
    boundaryHe->face = f;
  }
  centerVert.ptr->halfedge = trailingHalfedges[0];


  isCanonicalFlag = false;
  return centerVert;
}

Vertex HalfedgeMesh::collapseEdge(Edge e) {

  if (e.isBoundary()) {
    return collapseEdgeAlongBoundary(e);
  }

  // TODO for now only valid on triangle meshes and away from the boundary.

  // === Gather some elements

  Halfedge* heA0;
  // If there's a single boundary vertex, be sure we keep it
  if (e.halfedge().twin().vertex().isBoundary() && !e.halfedge().vertex().isBoundary()) {
    heA0 = e.halfedge().twin().ptr;
  } else {
    heA0 = e.halfedge().ptr;
  }

  Halfedge* heA1 = heA0->next;
  Halfedge* heA2 = heA1->next;
  Face* fA = heA0->face;
  Edge* eADelete = heA1->edge;
  Edge* eAKeep = heA2->edge;

  Halfedge* heB0 = heA0->twin;
  Halfedge* heB1 = heB0->next;
  Halfedge* heB2 = heB1->next;
  Face* fB = heB0->face;
  Edge* eBDelete = heB2->edge;
  Edge* eBKeep = heB1->edge;

  Vertex* vKeep = heA0->vertex;
  Vertex* vDiscard = heA0->twin->vertex;

  // === Check validity

  // collapsing around a degree-2 vertex can be done, but this code does not handle that correctly
  if (Vertex(vKeep).degree() <= 2 || Vertex(vDiscard).degree() <= 2) {
    return Vertex();
  }

  bool vKeepIsBoundary = Vertex(vKeep).isBoundary();

  // (should be exactly two vertices, the opposite diamond vertices, in the intersection of the 1-rings)
  std::unordered_set<Vertex> vKeepNeighbors;
  for (Vertex vN : Vertex(vKeep).adjacentVertices()) {
    vKeepNeighbors.insert(vN);
  }
  size_t nShared = 0;
  for (Vertex vN : Vertex(vDiscard).adjacentVertices()) {
    if (vKeepNeighbors.find(vN) != vKeepNeighbors.end()) {
      nShared++;
    }
  }
  if (nShared > 2) {
    return Vertex();
  }


  // === Update a whole bunch of pointers

  { // Update all of the halfedges around vDiscard (do this loop before we break things
    Halfedge* currHe = heA1;
    Halfedge* firstHe = heA1;
    do {
      currHe->vertex = vKeep;
      currHe = currHe->twin->next;
    } while (currHe != firstHe);
  }

  // Fix vertices
  if (vKeep->halfedge == heA0 || vKeep->halfedge == heB1) {
    // Only fix if needed, which ensure we don't mess up boundary vertices
    vKeep->halfedge = heB2->twin;
  }
  if (heA2->vertex->halfedge == heA2) {
    heA2->vertex->halfedge = heA1->twin;
  }
  if (heB2->vertex->halfedge == heB2) {
    heB2->vertex->halfedge = heB2->twin->next;
  }

  // Fix edges
  eAKeep->halfedge = heA2->twin;
  eBKeep->halfedge = heB1->twin;

  // Fix halfedges
  heA1->twin->edge = eAKeep;
  heA1->twin->twin = heA2->twin;
  heA2->twin->twin = heA1->twin;
  heB2->twin->edge = eBKeep;
  heB2->twin->twin = heB1->twin;
  heB1->twin->twin = heB2->twin;

  if (vKeepIsBoundary) {
    heA1->edge->isBoundary = true;
    heB2->edge->isBoundary = true;
  }


  // === Delete everything which needs to be deleted
  deleteElement(e);
  deleteElement(heA0);
  deleteElement(heA1);
  deleteElement(heA2);
  deleteElement(heB0);
  deleteElement(heB1);
  deleteElement(heB2);
  deleteElement(fA);
  deleteElement(fB);
  deleteElement(vDiscard);
  deleteElement(eADelete);
  deleteElement(eBDelete);

  isCanonicalFlag = false;

  return Vertex(vKeep);
}


Vertex HalfedgeMesh::collapseEdgeAlongBoundary(Edge e) {

  if (!e.isBoundary()) throw std::runtime_error("Called away from boundary");


  if (e.halfedge().next().edge().isBoundary() && e.halfedge().next().next().edge().isBoundary()) {
    throw std::runtime_error("Tried to collapse single face");
  }

  Halfedge* heA0 = e.halfedge().ptr;
  Halfedge* heA1 = heA0->next;
  Halfedge* heA2 = heA1->next;
  Face* fA = heA0->face;
  Edge* eADelete = heA1->edge;
  Edge* eAKeep = heA2->edge;

  Halfedge* heB = heA0->twin;
  Halfedge* heBNext = heB->next;
  Halfedge* heBPrev = heA0; // about to change
  while (heBPrev->next != heB) heBPrev = heBPrev->next->twin;
  Face* bL = heB->face;


  Vertex* vKeep = heA0->vertex;
  Vertex* vDiscard = heA0->twin->vertex;

  // (should be exactly two vertices, the opposite diamond vertices, in the intersection of the 1-rings)
  std::unordered_set<Vertex> vKeepNeighbors;
  for (Vertex vN : Vertex(vKeep).adjacentVertices()) {
    vKeepNeighbors.insert(vN);
  }
  size_t nShared = 0;
  for (Vertex vN : Vertex(vDiscard).adjacentVertices()) {
    if (vKeepNeighbors.find(vN) != vKeepNeighbors.end()) {
      nShared++;
    }
  }
  if (nShared > 2) {
    return Vertex();
    cout << "can't collapse: vertex neighborhoods are not distinct" << endl;
  }

  // === Update pointers

  { // Update all of the halfedges around vDiscard (do this loop before we break things
    Halfedge* currHe = heA1;
    Halfedge* firstHe = heA1;
    do {
      currHe->vertex = vKeep;
      currHe = currHe->twin->next;
    } while (currHe != firstHe);
  }


  if (heA2->vertex->halfedge == heA2) {
    heA2->vertex->halfedge = heA1->twin;
  }

  vKeep->halfedge = heBPrev->twin;

  // Fix edges
  if (heA2->twin->isReal) {
    eAKeep->halfedge = heA2->twin;
  } else {
    eAKeep->halfedge = heA1->twin;
  }

  // Fix halfedges
  heA1->twin->edge = eAKeep;
  heA1->twin->twin = heA2->twin;
  heA2->twin->twin = heA1->twin;
  heBPrev->next = heBNext;

  // Fix boundary loop
  bL->halfedge = heBPrev;

  ensureVertexHasBoundaryHalfedge(heA2->vertex);

  // === Delete everything which needs to be deleted
  deleteElement(e);
  deleteElement(heA0);
  deleteElement(heA1);
  deleteElement(heA2);
  deleteElement(heB);
  deleteElement(fA);
  deleteElement(vDiscard);
  deleteElement(eADelete);

  isCanonicalFlag = false;

  return Vertex(vKeep);
}

void HalfedgeMesh::ensureVertexHasBoundaryHalfedge(Vertex v) {
  if (!v.isBoundary()) return;
  Vertex* rawV = v.ptr;
  while (rawV->halfedge->twin->isReal) {
    rawV->halfedge = rawV->halfedge->twin->next;
  }
}

bool HalfedgeMesh::removeFaceAlongBoundary(Face f) {

  // Find the boundary halfedge
  Halfedge heBoundary;
  int bCount = 0;
  for (Halfedge he : f.adjacentHalfedges()) {
    if (!he.twin().isReal()) {
      bCount++;
      heBoundary = he;
    }
  }
  if (bCount == 0) {
    throw std::runtime_error("called on non-boundary face");
  }
  if (bCount == 1) {
    // Remove a non-ear boundary face with one boundary edge

    Halfedge* he0 = heBoundary.ptr;
    Halfedge* he0T = he0->twin;
    Halfedge* he1 = he0->next;
    Halfedge* he2 = he1->next;
    Vertex* v0 = he0->vertex;
    Vertex* v1 = he1->vertex;
    Vertex* v2 = he2->vertex;
    Face* fRemove = he0->face;
    Face* bLoop = he0T->face;

    // Vertex halfedges
    v0->halfedge = he2->twin;
    v2->halfedge = he1->twin;

    // Nexts
    he2->next = he0T->next;
    v1->halfedge->twin->next = he1;

    // Faces
    he1->face = bLoop;
    he2->face = bLoop;

    // mark boundary
    v2->isBoundary = true;
    he1->isReal = false;
    he2->isReal = false;

    deleteElement(he0->edge);
    deleteElement(he0);
    deleteElement(he0T);
    deleteElement(fRemove);

    isCanonicalFlag = false;
    return true;

  } else if (bCount == 2) {
    // Remove an "ear" along the boundary

    // Gather elements
    Halfedge* he0 = f.halfedge().ptr;
    while (!he0->twin->isReal) he0 = he0->next;
    Halfedge* he0T = he0->twin;
    Halfedge* he1 = he0->next;
    Halfedge* he1T = he1->twin;
    Edge* e1 = he1->edge;
    Halfedge* he2 = he1->next;
    Halfedge* he2T = he2->twin;
    Edge* e2 = he2->edge;
    Vertex* v0 = he0->vertex;
    Vertex* v1 = he1->vertex;
    Vertex* v2 = he2->vertex;
    Face* fRemove = he0->face;

    Halfedge* heNext = he1T->next;
    Halfedge* hePrev = he0T;
    while (hePrev->isReal) hePrev = hePrev->next->twin;

    // Vertex halfedges
    v0->halfedge = hePrev->twin;
    v1->halfedge = he0T;

    // Nexts
    hePrev->next = heNext;

    // Boundary loop
    hePrev->face->halfedge = hePrev;

    // mark boundary
    he0->isReal = false;

    deleteElement(fRemove);
    deleteElement(v2);
    deleteElement(he1);
    deleteElement(he1T);
    deleteElement(e1);
    deleteElement(he2);
    deleteElement(he2T);
    deleteElement(e2);

    isCanonicalFlag = false;
    return true;

  } else {
    // Remove entire component

    /*
    Halfedge* he0 = heBoundary.ptr;
    Halfedge* he0T = he0->twin;
    Edge* e0 = he0->edge;
    Halfedge* he1 = he0->next;
    Halfedge* he1T = he1->twin;
    Edge* e1 = he1->edge;
    Halfedge* he2 = he1->next;
    Halfedge* he2T = he2->twin;
    Edge* e2 = he2->edge;
    Vertex* v0 = he0->vertex;
    Vertex* v1 = he1->vertex;
    Vertex* v2 = he2->vertex;
    Face* fFace = he0->face;
    Face* fBound = he0T->face;


    deleteElement(he0);
    deleteElement(he1);
    deleteElement(he2);

    deleteElement(he0T);
    deleteElement(he1T);
    deleteElement(he2T);

    deleteElement(e0);
    deleteElement(e1);
    deleteElement(e2);

    deleteElement(v0);
    deleteElement(v1);
    deleteElement(v2);

    deleteElement(fFace);
    deleteElement(fBound);

    isCanonicalFlag = false;
    return true;
    */

    // The removal/insertion code doesn't support changing boundary structure yet
    return false;
  }
}

void HalfedgeMesh::setEdgeHalfedge(Edge e, Halfedge he) { e.ptr->halfedge = he.ptr; }

std::vector<Face> HalfedgeMesh::triangulate(Face f) {

  if (f.degree() == 3) {
    return {f};
  }


  std::vector<DynamicVertex> neighVerts;
  std::vector<DynamicFace> allFaces;

  for (Vertex v : f.adjacentVertices()) {
    neighVerts.emplace_back(v, this);
  }

  allFaces.emplace_back(f, this);
  Halfedge currHe = f.halfedge();
  for (size_t i = 2; i + 1 < neighVerts.size(); i++) {
    Halfedge newHe = connectVertices(currHe.face(), neighVerts[0], neighVerts[i]);
    allFaces.emplace_back(newHe.twin().face(), this);
    currHe = newHe;
  }

  std::vector<Face> staticFaces;
  for (DynamicFace& d : allFaces) {
    staticFaces.emplace_back(d);
  }

  isCanonicalFlag = false;
  return staticFaces;
}

void HalfedgeMesh::validateConnectivity() {

  // == Halfedges

  // Check valid pointers
  for (Halfedge& he : rawHalfedges) {
    if (he.isDead()) continue;
    if (he.twin == nullptr || !&*he.twin || he.twin->isDead()) throw std::logic_error("bad twin pointer");
    if (he.next == nullptr || !&*he.next || he.next->isDead()) throw std::logic_error("bad next pointer");
    if (he.vertex == nullptr || !&*he.vertex || he.vertex->isDead()) throw std::logic_error("bad vertex pointer");
    if (he.edge == nullptr || !&*he.edge || he.edge->isDead()) throw std::logic_error("bad edge pointer");
    if (he.face == nullptr || !&*he.face || he.face->isDead()) throw std::logic_error("bad face pointer");
  }
  for (Vertex& v : rawVertices) {
    if (v.isDead()) continue;
    if (v.halfedge == nullptr || !&*v.halfedge || v.halfedge->isDead())
      throw std::logic_error("bad halfedge pointer in vertex");
  }
  for (Edge& e : rawEdges) {
    if (e.isDead()) continue;
    if (e.halfedge == nullptr || !&*e.halfedge || e.halfedge->isDead())
      throw std::logic_error("bad halfedge pointer in edge");
  }
  for (Face& f : rawFaces) {
    if (f.isDead()) continue;
    if (f.halfedge == nullptr || !&*f.halfedge || f.halfedge->isDead())
      throw std::logic_error("bad halfedge pointer in face");
  }

  // Check edge and twin sanity
  for (Halfedge& he : rawHalfedges) {
    if (he.isDead()) continue;
    if (&he != he.twin->twin) throw std::logic_error("twins not reflective");
    if (&he == he.twin) throw std::logic_error("self-twin");
    if (&he != he.edge->halfedge && he.twin != he.edge->halfedge)
      throw std::logic_error("edge.halfedge doesn't match halfedge.edge");
  }

  // Check face & next sanity
  for (Face& f : rawFaces) {
    if (f.isDead()) continue;
    Halfedge* currHe = f.halfedge;
    Halfedge* firstHe = f.halfedge;
    size_t count = 0;
    do {
      if (currHe->face != &f) throw std::logic_error("face.halfedge doesn't match halfedge.face");
      currHe = currHe->next;
      count++;
      if (count > rawHalfedges.size()) throw std::logic_error("next forms non-face loop");
    } while (currHe != firstHe);

    if (count < 3) throw std::logic_error("face of degree < 2");
  }

  for (Face& f : rawBoundaryLoops) {
    if (f.isDead()) continue;
    Halfedge* currHe = f.halfedge;
    Halfedge* firstHe = f.halfedge;
    size_t count = 0;
    do {
      if (currHe->face != &f) throw std::logic_error("(boundary loop) face.halfedge doesn't match halfedge.face");
      currHe = currHe->next;
      count++;
      if (count > rawHalfedges.size()) throw std::logic_error("(boundary loop) next forms non-face loop");
    } while (currHe != firstHe);

    if (count < 3) throw std::logic_error("(boundary loop) face of degree < 2");
  }

  for (Halfedge& he : rawHalfedges) {
    if (he.isDead()) continue;
    Halfedge* currHe = &he;
    Halfedge* firstHe = &he;
    size_t count = 0;
    do {
      if (currHe->face != he.face) throw std::logic_error("he.next.**.face doesn't match he.face");
      currHe = currHe->next;
      count++;
      if (count > rawHalfedges.size()) throw std::logic_error("next forms non-face loop");
    } while (currHe != firstHe);

    // This can happen in irregular triangulations
    // if (he.vertex == he.next->twin->vertex) throw std::logic_error("halfedge face spur");

    if (he.vertex != he.twin->next->vertex) throw std::logic_error("halfedge vertices don't match");
  }

  // Check vertex orbit sanity
  for (Vertex& v : rawVertices) {
    if (v.isDead()) continue;
    Halfedge* currHe = v.halfedge;
    Halfedge* firstHe = v.halfedge;
    size_t count = 0;
    do {
      if (currHe->vertex != &v) throw std::logic_error("vertex.halfedge doesn't match halfedge.vertex");
      currHe = currHe->twin->next;
      count++;
      if (count > rawHalfedges.size()) throw std::logic_error("twin->next forms non-vertex loop");
    } while (currHe != firstHe);
  }

  // Verify isBoundary flag is correct
  for (Vertex& v : rawVertices) {
    if (v.isDead()) continue;
    bool hasBoundaryHe = false;
    Halfedge* currHe = v.halfedge;
    Halfedge* firstHe = v.halfedge;
    do {
      if (!currHe->isReal) hasBoundaryHe = true;
      currHe = currHe->twin->next;
    } while (currHe != firstHe);

    if (hasBoundaryHe != v.isBoundary) {
      throw std::logic_error("cached v.isBoundary is wrong");
    }
  }

  // Verify that for a mesh with boundary, outgoing boundary halfedge is the one that starts the boundary wedge
  // (aka it's the unique real halfedge along the boundary)
  for (Vertex& v : rawVertices) {
    if (v.isDead() || !v.isBoundary) continue;

    if (!v.halfedge->isReal) throw std::logic_error("v.halfedge() is not real");
    if (v.halfedge->twin->isReal) throw std::logic_error("v.halfedge() is not along boundary");
  }
}


Halfedge* HalfedgeMesh::getNewHalfedge(bool real) {

  // The boring case, when no resize is needed
  if (rawHalfedges.size() < rawHalfedges.capacity()) {
    rawHalfedges.emplace_back();
  }
  // The intesting case, where the vector resizes and we need to update pointers.
  else {

    Halfedge* oldStart = &rawHalfedges.front();

    // === Prep the "before" lists
    std::vector<std::ptrdiff_t> offsetsTwin(rawHalfedges.size());
    std::vector<std::ptrdiff_t> offsetsNext(rawHalfedges.size());
    std::vector<std::ptrdiff_t> offsetsV(rawVertices.size());
    std::vector<std::ptrdiff_t> offsetsE(rawEdges.size());
    std::vector<std::ptrdiff_t> offsetsF(rawFaces.size());
    std::vector<std::ptrdiff_t> offsetsB(rawBoundaryLoops.size());

    for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
      if (!rawHalfedges[iHe].isDead()) {
        offsetsTwin[iHe] = rawHalfedges[iHe].twin - oldStart;
        offsetsNext[iHe] = rawHalfedges[iHe].next - oldStart;
      }
    }
    for (size_t iV = 0; iV < rawVertices.size(); iV++) {
      if (!rawVertices[iV].isDead()) {
        offsetsV[iV] = rawVertices[iV].halfedge - oldStart;
      }
    }
    for (size_t iE = 0; iE < rawEdges.size(); iE++) {
      if (!rawEdges[iE].isDead()) {
        offsetsE[iE] = rawEdges[iE].halfedge - oldStart;
      }
    }
    for (size_t iF = 0; iF < rawFaces.size(); iF++) {
      if (!rawFaces[iF].isDead()) {
        offsetsF[iF] = rawFaces[iF].halfedge - oldStart;
      }
    }
    for (size_t iB = 0; iB < rawBoundaryLoops.size(); iB++) {
      if (!rawBoundaryLoops[iB].isDead()) {
        offsetsB[iB] = rawBoundaryLoops[iB].halfedge - oldStart;
      }
    }


    // Create a new halfedge, allowing the list to expand
    rawHalfedges.emplace_back();
    Halfedge* newStart = &rawHalfedges.front();

    // === Loop back through, shifting all pointers
    for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
      if (!rawHalfedges[iHe].isDead()) {
        rawHalfedges[iHe].twin = newStart + offsetsTwin[iHe];
        rawHalfedges[iHe].next = newStart + offsetsNext[iHe];
      }
    }
    for (size_t iV = 0; iV < rawVertices.size(); iV++) {
      if (!rawVertices[iV].isDead()) {
        rawVertices[iV].halfedge = newStart + offsetsV[iV];
      }
    }
    for (size_t iE = 0; iE < rawEdges.size(); iE++) {
      if (!rawEdges[iE].isDead()) {
        rawEdges[iE].halfedge = newStart + offsetsE[iE];
      }
    }
    for (size_t iF = 0; iF < rawFaces.size(); iF++) {
      if (!rawFaces[iF].isDead()) {
        rawFaces[iF].halfedge = newStart + offsetsF[iF];
      }
    }
    for (size_t iB = 0; iB < rawBoundaryLoops.size(); iB++) {
      if (!rawBoundaryLoops[iB].isDead()) {
        rawBoundaryLoops[iB].halfedge = newStart + offsetsB[iB];
      }
    }

    // Invoke relevant callback functions
    for (auto& f : halfedgeExpandCallbackList) {
      f(rawHalfedges.capacity());
    }
  }

  rawHalfedges.back().ID = nextElemID++;
  rawHalfedges.back().isReal = real;
  rawHalfedges.back().markDead(); // temporarily, to ensure we don't follow pointers
  if (real) {
    nRealHalfedgesCount++;
  } else {
    nImaginaryHalfedgesCount++;
  }
#ifndef NDEBUG
  rawHalfedges.back().parentMesh = this;
#endif
  return &rawHalfedges.back();
}

Vertex* HalfedgeMesh::getNewVertex() {

  // The boring case, when no resize is needed
  if (rawVertices.size() < rawVertices.capacity()) {
    rawVertices.emplace_back();
  }
  // The intesting case, where the vector resizes and we need to update pointers.
  else {

    Vertex* oldStart = &rawVertices.front();
    std::vector<std::ptrdiff_t> offsets(rawHalfedges.size());
    for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
      if (!rawHalfedges[iHe].isDead()) {
        offsets[iHe] = rawHalfedges[iHe].vertex - oldStart;
      }
    }

    // Create a new element, allowing the list to expand
    rawVertices.emplace_back();

    Vertex* newStart = &rawVertices.front();
    for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
      if (!rawHalfedges[iHe].isDead()) {
        rawHalfedges[iHe].vertex = newStart + offsets[iHe];
      }
    }

    // Invoke relevant callback functions
    for (auto& f : vertexExpandCallbackList) {
      f(rawVertices.capacity());
    }
  }

  rawVertices.back().ID = nextElemID++;
  rawVertices.back().markDead(); // temporarily, to ensure we don't follow pointers
  nVerticesCount++;
#ifndef NDEBUG
  rawVertices.back().parentMesh = this;
#endif
  return &rawVertices.back();
}

Edge* HalfedgeMesh::getNewEdge() {

  // The boring case, when no resize is needed
  if (rawEdges.size() < rawEdges.capacity()) {
    rawEdges.emplace_back();
  }
  // The intesting case, where the vector resizes and we need to update pointers.
  else {

    Edge* oldStart = &rawEdges.front();
    std::vector<std::ptrdiff_t> offsets(rawHalfedges.size());
    for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
      if (!rawHalfedges[iHe].isDead()) {
        offsets[iHe] = rawHalfedges[iHe].edge - oldStart;
      }
    }

    // Create a new element, allowing the list to expand
    rawEdges.emplace_back();

    Edge* newStart = &rawEdges.front();
    for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
      if (!rawHalfedges[iHe].isDead()) {
        rawHalfedges[iHe].edge = newStart + offsets[iHe];
      }
    }

    // Invoke relevant callback functions
    for (auto& f : edgeExpandCallbackList) {
      f(rawEdges.capacity());
    }
  }

  rawEdges.back().ID = nextElemID++;
  rawEdges.back().markDead(); // temporarily, to ensure we don't follow pointers
  nEdgesCount++;
#ifndef NDEBUG
  rawEdges.back().parentMesh = this;
#endif
  return &rawEdges.back();
}

Face* HalfedgeMesh::getNewFace() {

  // The boring case, when no resize is needed
  if (rawFaces.size() < rawFaces.capacity()) {
    rawFaces.emplace_back();
  }
  // The intesting case, where the vector resizes and we need to update pointers.
  else {
    Face* oldStart = &rawFaces.front();
    std::vector<std::ptrdiff_t> offsets(rawHalfedges.size());
    for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
      if (!rawHalfedges[iHe].isDead()) {
        offsets[iHe] = rawHalfedges[iHe].face - oldStart;
      }
    }

    // Create new element, allowing the list to expand
    rawFaces.emplace_back();

    Face* newStart = &rawFaces.front();
    for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
      if (!rawHalfedges[iHe].isDead()) {
        rawHalfedges[iHe].face = newStart + offsets[iHe];
      }
    }

    // Invoke relevant callback functions
    for (auto& f : faceExpandCallbackList) {
      f(rawFaces.capacity());
    }
  }

  rawFaces.back().ID = nextElemID++;
  rawFaces.back().isReal = true;
  rawFaces.back().markDead(); // temporarily, to ensure we don't follow pointers
  nFacesCount++;
#ifndef NDEBUG
  rawFaces.back().parentMesh = this;
#endif
  return &rawFaces.back();
}

void HalfedgeMesh::deleteElement(Halfedge he) {
  he->markDead();
  isCompressedFlag = false;

  if (he.isReal()) {
    nRealHalfedgesCount--;
  } else {
    nImaginaryHalfedgesCount--;
  }
}

void HalfedgeMesh::deleteElement(Edge e) {
  e->markDead();
  isCompressedFlag = false;
  nEdgesCount--;
}

void HalfedgeMesh::deleteElement(Vertex v) {
  v->markDead();
  isCompressedFlag = false;
  nVerticesCount--;
}

void HalfedgeMesh::deleteElement(Face f) {
  f->markDead();
  isCompressedFlag = false;
  nFacesCount--;
}

void HalfedgeMesh::compressHalfedges() {

  // Build the compressing shift
  std::vector<size_t> newIndMap; // maps new ind -> old ind
  std::vector<size_t> oldIndMap; // maps new old -> new ind
  oldIndMap.resize(rawHalfedges.size());
  for (size_t i = 0; i < rawHalfedges.size(); i++) {
    if (!rawHalfedges[i].isDead()) {
      oldIndMap[i] = newIndMap.size();
      newIndMap.push_back(i);
    }
  }
  // === Prep the "before" lists
  Halfedge* oldStart = &rawHalfedges.front();
  std::vector<size_t> offsetsTwin(rawHalfedges.size());
  std::vector<size_t> offsetsNext(rawHalfedges.size());
  std::vector<size_t> offsetsV(rawVertices.size());
  std::vector<size_t> offsetsE(rawEdges.size());
  std::vector<size_t> offsetsF(rawFaces.size());
  std::vector<size_t> offsetsB(rawBoundaryLoops.size());

  for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
    if (!rawHalfedges[iHe].isDead()) {
      offsetsTwin[iHe] = rawHalfedges[iHe].twin - oldStart;
      offsetsNext[iHe] = rawHalfedges[iHe].next - oldStart;
    }
  }
  for (size_t iV = 0; iV < rawVertices.size(); iV++) {
    if (!rawVertices[iV].isDead()) {
      offsetsV[iV] = rawVertices[iV].halfedge - oldStart;
    }
  }
  for (size_t iE = 0; iE < rawEdges.size(); iE++) {
    if (!rawEdges[iE].isDead()) {
      offsetsE[iE] = rawEdges[iE].halfedge - oldStart;
    }
  }
  for (size_t iF = 0; iF < rawFaces.size(); iF++) {
    if (!rawFaces[iF].isDead()) {
      offsetsF[iF] = rawFaces[iF].halfedge - oldStart;
    }
  }
  for (size_t iB = 0; iB < rawBoundaryLoops.size(); iB++) {
    if (!rawBoundaryLoops[iB].isDead()) {
      offsetsB[iB] = rawBoundaryLoops[iB].halfedge - oldStart;
    }
  }

  // Apply the permutation
  rawHalfedges = applyPermutation(rawHalfedges, newIndMap);
  Halfedge* newStart = &rawHalfedges.front();

  // === Loop back through, shifting all pointers
  // TODO since this is in compress(), should never be dead, right?
  for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
    if (!rawHalfedges[iHe].isDead()) {
      rawHalfedges[iHe].twin = newStart + oldIndMap[offsetsTwin[newIndMap[iHe]]];
      rawHalfedges[iHe].next = newStart + oldIndMap[offsetsNext[newIndMap[iHe]]];
    }
  }
  for (size_t iV = 0; iV < rawVertices.size(); iV++) {
    if (!rawVertices[iV].isDead()) {
      rawVertices[iV].halfedge = newStart + oldIndMap[offsetsV[iV]];
    }
  }
  for (size_t iE = 0; iE < rawEdges.size(); iE++) {
    if (!rawEdges[iE].isDead()) {
      rawEdges[iE].halfedge = newStart + oldIndMap[offsetsE[iE]];
    }
  }
  for (size_t iF = 0; iF < rawFaces.size(); iF++) {
    if (!rawFaces[iF].isDead()) {
      rawFaces[iF].halfedge = newStart + oldIndMap[offsetsF[iF]];
    }
  }
  for (size_t iB = 0; iB < rawBoundaryLoops.size(); iB++) {
    if (!rawBoundaryLoops[iB].isDead()) {
      rawBoundaryLoops[iB].halfedge = newStart + oldIndMap[offsetsB[iB]];
    }
  }

  // Invoke callbacks
  for (auto& f : halfedgePermuteCallbackList) {
    f(newIndMap);
  }
}

void HalfedgeMesh::compressEdges() {

  // Build the compressing shift
  std::vector<size_t> newIndMap; // maps new ind -> old ind
  std::vector<size_t> oldIndMap; // maps new old -> new ind
  oldIndMap.resize(rawEdges.size());
  for (size_t i = 0; i < rawEdges.size(); i++) {
    if (!rawEdges[i].isDead()) {
      oldIndMap[i] = newIndMap.size();
      newIndMap.push_back(i);
    }
  }


  Edge* oldStart = &rawEdges.front();
  std::vector<size_t> offsets(rawHalfedges.size());
  for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
    if (!rawHalfedges[iHe].isDead()) {
      offsets[iHe] = rawHalfedges[iHe].edge - oldStart;
    }
  }

  // Apply the permutation
  rawEdges = applyPermutation(rawEdges, newIndMap);

  Edge* newStart = &rawEdges.front();
  for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
    if (!rawHalfedges[iHe].isDead()) {
      rawHalfedges[iHe].edge = newStart + oldIndMap[offsets[iHe]];
    }
  }

  // Invoke callbacks
  for (auto& f : edgePermuteCallbackList) {
    f(newIndMap);
  }
}

void HalfedgeMesh::compressFaces() {

  // Build the compressing shift
  std::vector<size_t> newIndMap; // maps new ind -> old ind
  std::vector<size_t> oldIndMap; // maps new old -> new ind
  oldIndMap.resize(rawFaces.size());
  for (size_t i = 0; i < rawFaces.size(); i++) {
    if (!rawFaces[i].isDead()) {
      oldIndMap[i] = newIndMap.size();
      newIndMap.push_back(i);
    }
  }

  Face* oldStart = &rawFaces.front();
  std::vector<size_t> offsets(rawHalfedges.size());
  for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
    if (!rawHalfedges[iHe].isDead()) {
      offsets[iHe] = rawHalfedges[iHe].face - oldStart;
    }
  }

  // Apply the permutation
  rawFaces = applyPermutation(rawFaces, newIndMap);

  Face* newStart = &rawFaces.front();
  for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
    if (!rawHalfedges[iHe].isDead() && rawHalfedges[iHe].isReal) {
      rawHalfedges[iHe].face = newStart + oldIndMap[offsets[iHe]];
    }
  }

  // Invoke callbacks
  for (auto& f : facePermuteCallbackList) {
    f(newIndMap);
  }
}


void HalfedgeMesh::compressVertices() {

  // Build the compressing shift
  std::vector<size_t> newIndMap; // maps new ind -> old ind
  std::vector<size_t> oldIndMap; // maps old ind -> new ind
  oldIndMap.resize(rawVertices.size());
  for (size_t i = 0; i < rawVertices.size(); i++) {
    if (!rawVertices[i].isDead()) {
      oldIndMap[i] = newIndMap.size();
      newIndMap.push_back(i);
    }
  }

  Vertex* oldStart = &rawVertices.front();
  std::vector<size_t> offsets(rawHalfedges.size());
  for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
    if (!rawHalfedges[iHe].isDead()) {
      offsets[iHe] = rawHalfedges[iHe].vertex - oldStart;
    }
  }

  // Apply the permutation
  rawVertices = applyPermutation(rawVertices, newIndMap);

  Vertex* newStart = &rawVertices.front();
  for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
    if (!rawHalfedges[iHe].isDead()) {
      rawHalfedges[iHe].vertex = newStart + oldIndMap[offsets[iHe]];
    }
  }

  // Invoke callbacks
  for (auto& f : vertexPermuteCallbackList) {
    f(newIndMap);
  }
}

void HalfedgeMesh::compress() {

  if (isCompressed()) {
    return;
  }

  compressHalfedges();
  compressEdges();
  compressFaces();
  compressVertices();
  isCompressedFlag = true;
}


void HalfedgeMesh::canonicalize() {

  if (isCanonical()) {
    return;
  }
  compress();

  // Vertices and faces are automatically canonical, nothing to do there.

  { // == Reorder halfedges

    std::vector<size_t> newIndMap(rawHalfedges.size(), 7777777777777); // maps new ind -> old ind
    std::vector<size_t> oldIndMap(rawHalfedges.size(), 7777777777777); // maps new old -> new ind

    // Compute the permutation to canonical
    Halfedge* start = &rawHalfedges.front();
    size_t nextInd = 0;
    for (Face f : faces()) {
      for (Halfedge he : f.adjacentHalfedges()) {
        size_t currInd = he.ptr - start;
        size_t newInd = nextInd;
        newIndMap[newInd] = currInd;
        oldIndMap[currInd] = newInd;
        nextInd++;
      }
    }
    for (BoundaryLoop b : boundaryLoops()) {
      for (Halfedge he : b.adjacentHalfedges()) {
        size_t currInd = he.ptr - start;
        size_t newInd = nextInd;
        newIndMap[newInd] = currInd;
        oldIndMap[currInd] = newInd;
        nextInd++;
      }
    }

    // Apply the permutation
    // (can't just assign, need to avoid moving memory)
    std::vector<Halfedge> pHalfedge = applyPermutation(rawHalfedges, newIndMap);
    for (size_t i = 0; i < pHalfedge.size(); i++) {
      rawHalfedges[i] = pHalfedge[i];
    }

    // Fix up pointers.
    for (Halfedge& he : rawHalfedges) {
      he.twin = start + oldIndMap[(he.twin - start)];
      he.next = start + oldIndMap[(he.next - start)];
    }
    for (Vertex& v : rawVertices) {
      v.halfedge = start + oldIndMap[(v.halfedge - start)];
    }
    for (Edge& e : rawEdges) {
      e.halfedge = start + oldIndMap[(e.halfedge - start)];
    }
    for (Face& f : rawFaces) {
      f.halfedge = start + oldIndMap[(f.halfedge - start)];
    }
    for (Face& f : rawBoundaryLoops) {
      f.halfedge = start + oldIndMap[(f.halfedge - start)];
    }

    // Invoke callbacks
    for (auto& f : halfedgePermuteCallbackList) {
      f(newIndMap);
    }
  }


  { // == Reorder edges

    std::vector<size_t> newIndMap(rawEdges.size(), 7777777777777); // maps new ind -> old ind
    std::vector<size_t> oldIndMap(rawEdges.size(), 7777777777777); // maps new old -> new ind
    EdgeData<char> edgeSeen(this, false);

    // Compute the permutation to canonical
    Edge* start = &rawEdges.front();
    size_t nextInd = 0;
    for (Face f : faces()) {
      for (Halfedge he : f.adjacentHalfedges()) {
        if (!edgeSeen[he.edge()]) {
          size_t currInd = he.edge().ptr - start;
          size_t newInd = nextInd;
          newIndMap[newInd] = currInd;
          oldIndMap[currInd] = newInd;
          nextInd++;
          edgeSeen[he.edge()] = true;
        }
      }
    }


    // Apply the permutation
    // (can't just assign, need to avoid moving memory)
    std::vector<Edge> pEdge = applyPermutation(rawEdges, newIndMap);
    for (size_t i = 0; i < pEdge.size(); i++) {
      rawEdges[i] = pEdge[i];
    }

    // Fix up pointers.
    for (Halfedge& he : rawHalfedges) {
      he.edge = start + oldIndMap[(he.edge - start)];
    }

    // Invoke callbacks
    for (auto& f : edgePermuteCallbackList) {
      f(newIndMap);
    }
  }

  isCanonicalFlag = true;
}

size_t HalfedgeMesh::indexOf(Halfedge* ptr) { return (ptr - &rawHalfedges[0]); }
size_t HalfedgeMesh::indexOf(Vertex* ptr) { return (ptr - &rawVertices[0]); }
size_t HalfedgeMesh::indexOf(Edge* ptr) { return (ptr - &rawEdges[0]); }
size_t HalfedgeMesh::indexOf(Face* ptr) { return (ptr - &rawFaces[0]); }


} // namespace halfedge_mesh
} // namespace geometrycentral
