#pragma once

#include <cmath>
#include <utility>
#include <queue>
#include <vector>

#include "geometrycentral/geometry/geometry.h"
#include "geometrycentral/geometry/edge_length_geometry.h"
#include "geometrycentral/utilities/utilities.h"

using namespace geometrycentral;

const double REL_ERR=1e-8;

struct Window {
  HalfedgePtr halfedge;
  double b0, b1;
  double d0, d1;

  double pseudoSrcDist, minDist;
  VertexPtr src, pseudoSrc;
  int pseudoSrcBirthTime;
  int level;

  bool operator<(const Window& right) const { return minDist > right.minDist; }
  void computeMinDist() {
    double wLen = b1-b0;
    double xProj = (d0*d0+wLen*wLen-d1*d1)/(2.*wLen);
    if( xProj < 0. ) minDist = d0 + pseudoSrcDist;
    else if( xProj > wLen ) minDist = d1 + pseudoSrcDist;
    else minDist = sqrt(fabs(d0*d0-xProj*xProj)) + pseudoSrcDist;
  }
  Vector2 flattenedSrc() const {
    Vector2 src2D;
    double wLen = b1-b0;
    src2D.x = (d0*d0+wLen*wLen-d1*d1)/(2.*wLen);
    src2D.y = sqrt(fabs(d0*d0-src2D.x*src2D.x));
    src2D.x += b0;
    return src2D;
  }
};

struct PseudoWindow {
  VertexPtr v;
  double dist;
  VertexPtr src, pseudoSrc;
  unsigned pseudoSrcBirthTime;
  unsigned level;

  bool operator<(const PseudoWindow& right) const { return dist > right.dist; }
};

struct SplitInfo {
  double dist;
  VertexPtr pseudoSrc, src;
  unsigned level;
  double x;

  SplitInfo() {
    dist = std::numeric_limits<double>::infinity();
    pseudoSrc = nullptr;
    src = nullptr;
    level = -1;
    x = std::numeric_limits<double>::infinity();
  }
};

struct VertInfo {
  int birthTime;
  double dist;
  bool isSource;
  HalfedgePtr enterHalfedge;
  VertexPtr pseudoSrc, src;

  VertInfo() {
    birthTime = -1;
    dist = std::numeric_limits<double>::infinity();
    isSource = false;
    enterHalfedge = nullptr;
    pseudoSrc = nullptr;
    src = nullptr;
  }
};

struct GeodesicKeyPoint {
  bool isVertex;
  unsigned id;
  double pos;
};

class ExactPolyhedralGeodesics {

public:
  ExactPolyhedralGeodesics( EdgeLengthGeometry* geom_);

  void addSource(VertexPtr v);
  VertexData<double> computeDistance(VertexPtr v);
  VertexData<double> computeDistance();

private:
  HalfedgeMesh* mesh;
  EdgeLengthGeometry* geom;

  std::vector<VertexPtr> srcVerts;
  HalfedgeData<SplitInfo> splitInfos;
  VertexData<VertInfo> vertInfos;
  std::priority_queue<Window> winQ;
  std::priority_queue<PseudoWindow> pseudoSrcQ;

  std::vector<Window> storedWindows;
  std::vector<FacePtr> keptFaces;
  HalfedgeData<std::list<Window>> allWindows;
  bool keptAllWindows = false;

  // stats
  int numOfWinGen;
  int maxWinQSize, maxPseudoQSize;
  int totalCalcVertNum;
  double geodesicRadius;
  bool geodesicRadiusReached;

  void clear();
  void initialize();
  void propogateWindow(const Window& win);
  void generateSubWinsForPseudoSrc(const PseudoWindow& pseudoWin);
  void generateSubWinsForPseudoSrcFromWindow(const PseudoWindow& pseudoWin, HalfedgePtr& startHe, HalfedgePtr& endHe);
  void generateSubWinsForPseudoSrcFromPseudoSrc(const PseudoWindow& pseudoWin, HalfedgePtr& startHe, HalfedgePtr& endHe);
  void buildWindow(const Window& pWin, HalfedgePtr& he, double t0, double t1, const Vector2& v0, const Vector2& v1, Window& win);
  bool isValidWindow(const Window& win, bool isLeftChild);
  double intersect(const Vector2& v0, const Vector2& v1, const Vector2& p0, const Vector2& p1);
};
