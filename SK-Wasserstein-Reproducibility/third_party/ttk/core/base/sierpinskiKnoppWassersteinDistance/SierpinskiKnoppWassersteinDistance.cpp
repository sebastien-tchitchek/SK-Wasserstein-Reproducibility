#include <SierpinskiKnoppWassersteinDistance.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>

namespace {

struct PtLD {
  long double x, y;
};

static inline PtLD sub(const PtLD &a, const PtLD &b) {
  return PtLD{a.x - b.x, a.y - b.y};
}

static inline long double cross2(const PtLD &u, const PtLD &v) {
  return u.x * v.y - u.y * v.x;
}

static inline PtLD midpoint(const PtLD &a, const PtLD &b) {
  return PtLD{(a.x + b.x) * 0.5L, (a.y + b.y) * 0.5L};
}

static bool pointInTriangle(const PtLD &p,
                            const PtLD &a,
                            const PtLD &b,
                            const PtLD &c,
                            const long double eps) {
  auto s = [&](const PtLD &p1, const PtLD &p2, const PtLD &p3) {
    return cross2(sub(p1, p3), sub(p2, p3));
  };

  const long double s1 = s(p, a, b);
  const long double s2 = s(p, b, c);
  const long double s3 = s(p, c, a);

  const bool has_neg = (s1 < -eps) || (s2 < -eps) || (s3 < -eps);
  const bool has_pos = (s1 > eps) || (s2 > eps) || (s3 > eps);
  return !(has_neg && has_pos);
}

struct FractionTri {
  PtLD start, end, mid;
};

static inline FractionTri baseTriangle() {
  return FractionTri{PtLD{0.0L, 0.0L}, PtLD{1.0L, 1.0L}, PtLD{0.0L, 1.0L}};
}

static inline std::pair<FractionTri, FractionTri> splitOnce(
  const FractionTri &T) {
  const PtLD k = midpoint(T.start, T.end);
  FractionTri T0{T.start, T.mid, k};
  FractionTri T1{T.mid, T.end, k};
  return {T0, T1};
}

static inline int clampLevel(int L) {
  if(L < 0)
    return 0;
  if(L > 50)
    return 50;
  return L;
}

static inline long double clamp01(long double v) {
  if(v < 0.0L)
    return 0.0L;
  if(v > 1.0L)
    return 1.0L;
  return v;
}

static std::uint64_t skXYToKPrefix(const PtLD &p,
                                   const int L,
                                   const long double eps) {
  FractionTri T = baseTriangle();
  std::uint64_t k = 0;

  for(int i = 0; i < L; ++i) {
    auto children = splitOnce(T);
    const FractionTri &T0 = children.first;
    const FractionTri &T1 = children.second;

    int bit = 0;
    if(pointInTriangle(p, T0.start, T0.end, T0.mid, eps)) {
      bit = 0;
      T = T0;
    } else {
      bit = 1;
      T = T1;
    }

    k = (k << 1) | static_cast<std::uint64_t>(bit);
  }

  return k;
}

static double skSelectorLeft(int L, long double x, long double y) {
  L = clampLevel(L);

  x = clamp01(x);
  y = clamp01(y);
  if(x > y)
    std::swap(x, y);

  const PtLD p{x, y};
  const std::uint64_t k = skXYToKPrefix(p, L, 1e-15L);
  const long double denom = std::ldexp(1.0L, L);
  return static_cast<double>(static_cast<long double>(k) / denom);
}

template <typename PairType>
static inline bool finitePairToNormalizedXY(const PairType &pair,
                                            long double &x,
                                            long double &y) {
  x = static_cast<long double>(pair.birth.sfValue);
  y = static_cast<long double>(pair.death.sfValue);

  if(!std::isfinite(static_cast<double>(x))
     || !std::isfinite(static_cast<double>(y))) {
    return false;
  }

  x = clamp01(x);
  y = clamp01(y);

  if(x > y)
    std::swap(x, y);

  return y > x;
}

} // namespace

ttk::SierpinskiKnoppWassersteinDistance::SierpinskiKnoppWassersteinDistance() {
  this->setDebugMsgPrefix("SierpinskiKnoppWassersteinDistance");
}

int ttk::SierpinskiKnoppWassersteinDistance::execute(
  const DiagramType &diagram0,
  const DiagramType &diagram1,
  std::vector<ProjectionPoint> &output,
  const int L) const {

  ttk::Timer tm{};
  output.clear();
  output.reserve(2 * (diagram0.size() + diagram1.size()));

  const auto addDiagram = [&](const DiagramType &diagram, const int diagramId) {
    for(size_t i = 0; i < diagram.size(); ++i) {
      long double x{}, y{};
      if(!finitePairToNormalizedXY(diagram[i], x, y)) {
        continue;
      }

      const double birth = static_cast<double>(x);
      const double death = static_cast<double>(y);
      const double persistence = death - birth;
      const double tOff = skSelectorLeft(L, x, y);

      ProjectionPoint off{};
      off.t = tOff;
      off.x = birth;
      off.y = death;
      off.birth = birth;
      off.death = death;
      off.persistence = persistence;
      off.diagramId = diagramId;
      off.isDiagonal = 0;
      off.pairId = static_cast<int>(i);
      off.category = 2 * diagramId;
      output.emplace_back(off);

      const long double m = 0.5L * (x + y);
      const double dm = static_cast<double>(m);
      const double tDiag = skSelectorLeft(L, m, m);

      ProjectionPoint diag{};
      diag.t = tDiag;
      diag.x = dm;
      diag.y = dm;
      diag.birth = birth;
      diag.death = death;
      diag.persistence = persistence;
      diag.diagramId = diagramId;
      diag.isDiagonal = 1;
      diag.pairId = static_cast<int>(i);
      diag.category = 2 * diagramId + 1;
      output.emplace_back(diag);
    }
  };

  addDiagram(diagram0, 0);
  addDiagram(diagram1, 1);

  std::sort(output.begin(), output.end(), [](const ProjectionPoint &a,
                                             const ProjectionPoint &b) {
    if(a.t < b.t)
      return true;
    if(a.t > b.t)
      return false;
    if(a.diagramId < b.diagramId)
      return true;
    if(a.diagramId > b.diagramId)
      return false;
    if(a.isDiagonal < b.isDiagonal)
      return true;
    if(a.isDiagonal > b.isDiagonal)
      return false;
    return a.pairId < b.pairId;
  });

  this->printMsg("Projected " + std::to_string(output.size())
                   + " atoms on the SK line",
                 1.0, tm.getElapsedTime(), this->threadNumber_);

  return 1;
}
