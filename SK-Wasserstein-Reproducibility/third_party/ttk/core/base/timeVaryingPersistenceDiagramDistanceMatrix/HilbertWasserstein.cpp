#include "HilbertWasserstein.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

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
                            long double eps) {
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

static inline std::pair<FractionTri, FractionTri> splitOnce(const FractionTri &T) {
  const PtLD K = midpoint(T.start, T.end);
  FractionTri T0{T.start, T.mid, K};
  FractionTri T1{T.mid, T.end, K};
  return {T0, T1};
}

static inline int clampL(int L) {
  if(L < 0) return 0;
  if(L > 50) return 50;
  return L;
}

static std::uint64_t sk_xy_to_kprefix(const PtLD &p, int L, long double eps) {
  FractionTri T = baseTriangle();
  std::uint64_t k = 0;

  for(int i = 0; i < L; i++) {
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

static inline long double clamp01_ld(long double v) {
  if(v < 0.0L) return 0.0L;
  if(v > 1.0L) return 1.0L;
  return v;
}

static inline double hilbert_tmid(int L, long double x, long double y) {
  L = clampL(L);

  x = clamp01_ld(x);
  y = clamp01_ld(y);
  if(x > y) std::swap(x, y);

  const PtLD p{x, y};
  const std::uint64_t k = sk_xy_to_kprefix(p, L, 1e-15L);
  const long double denom = std::ldexp(1.0L, L); // 2^L
  const long double tmid = (static_cast<long double>(k) + 0.5L) / denom;
  return static_cast<double>(tmid);
}

static inline double hilbert_tleft(int L, long double x, long double y) {
  L = clampL(L);

  x = clamp01_ld(x);
  y = clamp01_ld(y);
  if(x > y)
    std::swap(x, y);

  const PtLD p{x, y};
  const std::uint64_t k = sk_xy_to_kprefix(p, L, 1e-15L);
  const long double denom = std::ldexp(1.0L, L); // 2^L
  const long double tleft = static_cast<long double>(k) / denom;
  return static_cast<double>(tleft);
}

static inline double sk_selector_value(int L, long double x, long double y) {
  return hilbert_tleft(L, x, y);
  // return hilbert_tmid(L, x, y);
}

static inline double next_from_two_sorted(const std::vector<double> &a,
                                          std::size_t &ia,
                                          const std::vector<double> &b,
                                          std::size_t &ib) {
  assert(ia < a.size() || ib < b.size());

  if(ib >= b.size())
    return a[ia++];
  if(ia >= a.size())
    return b[ib++];
  if(a[ia] <= b[ib])
    return a[ia++];
  return b[ib++];
}

static inline double clamp01(double x) {
  if(x < 0.0) return 0.0;
  if(x > 1.0) return 1.0;
  return x;
}
static inline double d_to_boundary(double x) {
  x = clamp01(x);
  return std::min(x, 1.0 - x);
}

static double pot_figalli_1d_fast_L1_sorted_cost(const std::vector<double> &xs,
                                                 const std::vector<double> &ys) {
  const int n = (int)xs.size();
  const int m = (int)ys.size();
  if(n == 0 && m == 0) return 0.0;

  if(n == 0) {
    double cost = 0.0;
    for(double y : ys) cost += d_to_boundary(y);
    return cost;
  }
  if(m == 0) {
    double cost = 0.0;
    for(double x : xs) cost += d_to_boundary(x);
    return cost;
  }

  std::vector<double> z;
  z.reserve((size_t)n + (size_t)m + 2);

  auto push_unique = [&](double v) {
    if(z.empty() || z.back() != v) z.push_back(v);
  };

  push_unique(0.0);
  int i = 0, j = 0;
  while(i < n || j < m) {
    double v;
    if(j >= m || (i < n && xs[i] <= ys[j])) v = xs[i++];
    else                                     v = ys[j++];
    push_unique(v);
  }
  push_unique(1.0);

  auto count_at = [](const std::vector<double> &arr, int &p, double x) -> int {
    int c = 0;
    while(p < (int)arr.size() && arr[p] == x) { ++c; ++p; }
    return c;
  };

  int p_mu = 0, p_nu = 0;
  int64_t H = 0;
  H += (int64_t)count_at(xs, p_mu, z[0]);
  H -= (int64_t)count_at(ys, p_nu, z[0]);

  std::vector<double> wByVal((size_t)(n + m + 1), 0.0);

  for(size_t kk = 0; kk + 1 < z.size(); ++kk) {
    const double cur = z[kk];
    const double nxt = z[kk + 1];
    const double w = nxt - cur;
    if(w > 0.0) {
      const int64_t val = -H;
      const int64_t idx = val + n; 
      if(0 <= idx && idx < (int64_t)wByVal.size())
        wByVal[(size_t)idx] += w;
    }
    H += (int64_t)count_at(xs, p_mu, nxt);
    H -= (int64_t)count_at(ys, p_nu, nxt);
  }

  const double half = 0.5;
  double cum = 0.0;
  int64_t alpha = 0;
  bool found = false;
  for(int64_t val = -n; val <= m; ++val) {
    cum += wByVal[(size_t)(val + n)];
    if(cum + 1e-15 >= half) { alpha = val; found = true; break; }
  }
  if(!found) alpha = (int64_t)m;

  struct NodeC { double x; int64_t cnt; };
  std::deque<NodeC> supply, demand;

  if(alpha > 0) supply.push_back(NodeC{0.0, alpha});
  else if(alpha < 0) demand.push_back(NodeC{0.0, -alpha});

  auto consume = [&](NodeC &sp, NodeC &dm, double &cost) {
    const int64_t delta = std::min(sp.cnt, dm.cnt);
    if(delta <= 0) return;
    cost += (double)delta * std::abs(sp.x - dm.x);
    sp.cnt -= delta;
    dm.cnt -= delta;
  };

  double cost = 0.0;
  i = 0; j = 0;
  while(i < n || j < m) {
    const double nx = (i < n ? xs[i] : 2.0);
    const double ny = (j < m ? ys[j] : 2.0);

    if(nx <= ny) supply.push_back(NodeC{xs[i++], 1});
    else         demand.push_back(NodeC{ys[j++], 1});

    while(!supply.empty() && !demand.empty()) {
      NodeC &sp = supply.front();
      NodeC &dm = demand.front();
      consume(sp, dm, cost);
      if(sp.cnt == 0) supply.pop_front();
      if(dm.cnt == 0) demand.pop_front();
    }
  }

  const double B1 = 1.0;
  while(!supply.empty()) {
    cost += (double)supply.front().cnt * std::abs(supply.front().x - B1);
    supply.pop_front();
  }
  while(!demand.empty()) {
    cost += (double)demand.front().cnt * std::abs(demand.front().x - B1);
    demand.pop_front();
  }

  return cost;
}

struct CantorGapKey {
  uint64_t prefix;
  uint16_t k;
  bool operator==(const CantorGapKey &o) const { return prefix == o.prefix && k == o.k; }
};

struct CantorGapKeyHash {
  size_t operator()(const CantorGapKey &g) const noexcept {
    uint64_t x = g.prefix ^ (uint64_t(g.k) * 0x9e3779b97f4a7c15ULL);
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return (size_t)x;
  }
};

struct CantorLoc {
  bool inS = true;
  CantorGapKey key{};
  long double a = 0.0L, b = 1.0L;
};

static inline double clamp01_d(double v) {
  if(v < 0.0) return 0.0;
  if(v > 1.0) return 1.0;
  return v;
}

static inline uint64_t to_dyadic_num(double x, int Pbits) {
  if(x <= 0.0) return 0ULL;
  if(x >= 1.0) return (Pbits >= 63 ? ~0ULL : (1ULL << Pbits));
  const long double denom = std::ldexp(1.0L, Pbits);
  long double num = llroundl((long double)x * denom);
  if(num < 0) num = 0;
  if(num > denom) num = denom;
  return (uint64_t)num;
}

static CantorLoc locate_cantor_S_or_gap(double x, int Pbits, int L4max) {
  CantorLoc out;
  if(Pbits < 1) Pbits = 1;
  if(Pbits > 60) Pbits = 60;

  const int L4cap = (Pbits + 1) / 2;
  const int L4 = std::max(1, std::min(L4max, L4cap));

  const uint64_t denom = (1ULL << Pbits);
  uint64_t num = to_dyadic_num(x, Pbits);
  if(num == 0 || num == denom) {
    out.inS = true;
    return out;
  }

  uint64_t rem = num;
  const uint64_t mask = denom - 1;
  uint64_t prefix = 0;

  for(int k = 1; k <= L4; ++k) {
    __uint128_t tmp = ((__uint128_t)rem) << 2;
    uint64_t digit = (uint64_t)(tmp >> Pbits);
    rem = (uint64_t)(tmp & mask);

    if(digit == 1 || digit == 2) {
      if(digit == 1 && rem == 0) {
        out.inS = true;
        return out;
      }
      out.inS = false;
      out.key.prefix = prefix;
      out.key.k = (uint16_t)k;

      const uint64_t leftNum  = (prefix << 2) + 1ULL;
      const uint64_t rightNum = (prefix << 2) + 3ULL;
      out.a = std::ldexp((long double)leftNum,  -2 * k);
      out.b = std::ldexp((long double)rightNum, -2 * k);
      return out;
    }

    prefix = (prefix << 2) | digit;
    if(rem == 0) {
      out.inS = true;
      return out;
    }
  }

  out.inS = true;
  return out;
}

static inline int clampInt(int v, int lo, int hi) {
  return (v < lo) ? lo : (v > hi) ? hi : v;
}

static inline long double pow4_neg_ld(int L) {
  return std::ldexp(1.0L, -2 * L);
}

static std::vector<double> build_generators_Gn(int n, int L_eps) {
  // n = profondeur de G : G0, G1, G2, ...
  // L_eps = parametre epsilon : eps = 4^{-L_eps}
  n     = clampInt(n, 0, 20);
  L_eps = clampInt(L_eps, 2, 26);

  const long double eps = pow4_neg_ld(L_eps);
  const long double a   = 0.25L - eps;
  const long double b   = 0.75L + eps;

  std::vector<long double> cur;
  cur.reserve((std::size_t)1 << (n + 2));
  cur.push_back(0.0L);
  cur.push_back(a);
  cur.push_back(b);
  cur.push_back(1.0L);

  auto isEndpoint = [](long double t) -> bool {
    return t == 0.0L || t == 1.0L;
  };

  for(int it = 0; it < n; ++it) {
    std::vector<long double> nxt = cur;
    nxt.reserve(cur.size() + 2 * (cur.size() > 2 ? (cur.size() - 2) : 0));

    for(long double t : cur) {
      if(isEndpoint(t))
        continue;

      nxt.push_back(0.25L * t);

      nxt.push_back(0.75L + 0.25L * t);
    }

    std::sort(nxt.begin(), nxt.end());
    nxt.erase(std::unique(nxt.begin(), nxt.end()), nxt.end());
    cur.swap(nxt);
  }

  std::vector<double> out;
  out.reserve(cur.size());
  for(long double t : cur)
    out.push_back((double)t);

  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

static double pot_figalli_1d_fast_L1_generatorsGn_cost_sorted(
  const std::vector<double> &t_sorted,
  const std::vector<double> &s_sorted,
  int GnLevel,
  int L_eps) {

  const std::vector<double> G = build_generators_Gn(GnLevel, L_eps);
  if(G.size() < 2)
    return 0.0;

  const std::size_t n = t_sorted.size();
  const std::size_t m = s_sorted.size();

  std::size_t i = 0, j = 0;
  double totalCost = 0.0;

  for(std::size_t g = 0; g + 1 < G.size(); ++g) {
    const double a = G[g];
    const double b = G[g + 1];
    const double len = b - a;
    if(!(len > 0.0))
      continue;

    while(i < n && t_sorted[i] <= a) ++i;
    const std::size_t i0 = i;
    while(i < n && t_sorted[i] <  b) ++i;
    const std::size_t i1 = i;

    while(j < m && s_sorted[j] <= a) ++j;
    const std::size_t j0 = j;
    while(j < m && s_sorted[j] <  b) ++j;
    const std::size_t j1 = j;

    if(i0 == i1 && j0 == j1)
      continue;

    std::vector<double> tLoc;
    tLoc.reserve(i1 - i0);
    for(std::size_t p = i0; p < i1; ++p) {
      tLoc.push_back((t_sorted[p] - a) / len);
    }

    std::vector<double> sLoc;
    sLoc.reserve(j1 - j0);
    for(std::size_t p = j0; p < j1; ++p) {
      sLoc.push_back((s_sorted[p] - a) / len);
    }

    const double cLoc = pot_figalli_1d_fast_L1_sorted_cost(tLoc, sLoc);
    totalCost += cLoc * len;
  }

  return totalCost;
}


static double pot_figalli_1d_fast_L1_generatorsS_cost_sorted(const std::vector<double> &t_sorted,
                                                             const std::vector<double> &s_sorted,
                                                             int Pbits,
                                                             int L4max) {
  const int n = (int)t_sorted.size();
  const int m = (int)s_sorted.size();
  if(n == 0 && m == 0) return 0.0;

  struct GapBucketCost {
    long double a = 0.0L, b = 1.0L;
    std::vector<double> mu_x;
    std::vector<double> nu_x;
  };

  std::unordered_map<CantorGapKey, GapBucketCost, CantorGapKeyHash> buckets;
  buckets.reserve((size_t)(n + m));

  auto add_mu = [&](double x) {
    x = clamp01_d(x);
    CantorLoc loc = locate_cantor_S_or_gap(x, Pbits, L4max);
    if(loc.inS) return;
    auto &B = buckets[loc.key];
    if(B.mu_x.empty() && B.nu_x.empty()) { B.a = loc.a; B.b = loc.b; }
    B.mu_x.push_back(x);
  };
  auto add_nu = [&](double x) {
    x = clamp01_d(x);
    CantorLoc loc = locate_cantor_S_or_gap(x, Pbits, L4max);
    if(loc.inS) return;
    auto &B = buckets[loc.key];
    if(B.mu_x.empty() && B.nu_x.empty()) { B.a = loc.a; B.b = loc.b; }
    B.nu_x.push_back(x);
  };

  for(double x : t_sorted) add_mu(x);
  for(double y : s_sorted) add_nu(y);

  long double totalCost = 0.0L;

  for(const auto &kv : buckets) {
    const GapBucketCost &B = kv.second;
    const long double a = B.a;
    const long double b = B.b;
    const long double len = (b - a);
    if(len <= 0.0L) continue;

    std::vector<double> tLoc, sLoc;
    tLoc.reserve(B.mu_x.size());
    sLoc.reserve(B.nu_x.size());

    for(double x : B.mu_x) {
      long double u = ((long double)x - a) / len;
      tLoc.push_back(clamp01_d((double)u));
    }
    for(double y : B.nu_x) {
      long double u = ((long double)y - a) / len;
      sLoc.push_back(clamp01_d((double)u));
    }

    const double cLoc = pot_figalli_1d_fast_L1_sorted_cost(tLoc, sLoc);
    totalCost += (long double)cLoc * len;
  }

  return (double)totalCost;
}

struct SkAugAtom {
  bool onDiagonal;
  double x;
  double y;
};

static inline bool sk_off_atom_less(const SkOffAtom &a, const SkOffAtom &b) {
  if(a.t < b.t)
    return true;
  if(a.t > b.t)
    return false;
  if(a.x < b.x)
    return true;
  if(a.x > b.x)
    return false;
  return a.y < b.y;
}

static inline SkAugAtom next_from_two_sorted_augmented(
  const std::vector<SkOffAtom> &off_sorted,
  std::size_t &i_off,
  const std::vector<double> &diag_sorted,
  std::size_t &i_diag) {

  assert(i_off < off_sorted.size() || i_diag < diag_sorted.size());

  if(i_diag >= diag_sorted.size()) {
    const auto &a = off_sorted[i_off++];
    return SkAugAtom{false, a.x, a.y};
  }

  if(i_off >= off_sorted.size()) {
    ++i_diag;
    return SkAugAtom{true, 0.0, 0.0};
  }

  if(off_sorted[i_off].t <= diag_sorted[i_diag]) {
    const auto &a = off_sorted[i_off++];
    return SkAugAtom{false, a.x, a.y};
  }

  ++i_diag;
  return SkAugAtom{true, 0.0, 0.0};
}

static inline double diagonal_distance_sq_from_off_atom(const SkAugAtom &a) {
  assert(!a.onDiagonal);
  const double pers = a.y - a.x;
  return 0.5 * pers * pers; // dist^2((x,y), Delta)
}

static inline double c_delta_sq(const SkAugAtom &a, const SkAugAtom &b) {
  if(a.onDiagonal && b.onDiagonal)
    return 0.0;

  if(a.onDiagonal && !b.onDiagonal)
    return diagonal_distance_sq_from_off_atom(b);

  if(!a.onDiagonal && b.onDiagonal)
    return diagonal_distance_sq_from_off_atom(a);

  const double dx = a.x - b.x;
  const double dy = a.y - b.y;
  return dx * dx + dy * dy;
}

}

void diagram_to_skot_sorted_vectors(const ttk::DiagramType &D,
                                    int L,
                                    std::vector<double> &out_off_sorted,
                                    std::vector<double> &out_diag_sorted) {
  out_off_sorted.clear();
  out_diag_sorted.clear();

  out_off_sorted.reserve(D.size());
  out_diag_sorted.reserve(D.size());

  auto clamp01_tt = [](long double v) -> long double {
    if(v < 0.0L)
      return 0.0L;
    if(v > 1.0L)
      return 1.0L;
    return v;
  };

  for(const auto &pt : D) {
    const long double b = static_cast<long double>(pt.birth.sfValue);
    const long double d = static_cast<long double>(pt.death.sfValue);
    if(!std::isfinite(static_cast<double>(b)) || !std::isfinite(static_cast<double>(d)))
      continue;

    long double x = clamp01_tt(b);
    long double y = clamp01_tt(d);
    if(x > y)
      std::swap(x, y);

    if(y <= x)
      continue;

    out_off_sorted.push_back(sk_selector_value(L, x, y));

    const long double m = 0.5L * (x + y);
    out_diag_sorted.push_back(sk_selector_value(L, m, m));
  }

  std::sort(out_off_sorted.begin(), out_off_sorted.end());
  std::sort(out_diag_sorted.begin(), out_diag_sorted.end());
}

double skot_distance_from_precomputed_sorted(
  const std::vector<double> &off1_sorted,
  const std::vector<double> &diag1_sorted,
  const std::vector<double> &off2_sorted,
  const std::vector<double> &diag2_sorted) {

  assert(off1_sorted.size() == diag1_sorted.size());
  assert(off2_sorted.size() == diag2_sorted.size());

  const std::size_t Nleft  = off1_sorted.size() + diag2_sorted.size();
  const std::size_t Nright = off2_sorted.size() + diag1_sorted.size();
  assert(Nleft == Nright);

  std::size_t i_off1 = 0, i_diag2 = 0;
  std::size_t i_off2 = 0, i_diag1 = 0;

  double cost = 0.0;
  for(std::size_t k = 0; k < Nleft; ++k) {
    const double a = next_from_two_sorted(off1_sorted, i_off1, diag2_sorted, i_diag2);
    const double b = next_from_two_sorted(off2_sorted, i_off2, diag1_sorted, i_diag1);
    cost += std::abs(a - b);
  }

  return cost;
}

double skot_distance_normalized(const ttk::DiagramType &D1,
                                const ttk::DiagramType &D2,
                                int L) {
  std::vector<double> off1_sorted, diag1_sorted;
  std::vector<double> off2_sorted, diag2_sorted;

  diagram_to_skot_sorted_vectors(D1, L, off1_sorted, diag1_sorted);
  diagram_to_skot_sorted_vectors(D2, L, off2_sorted, diag2_sorted);

  return skot_distance_from_precomputed_sorted(
    off1_sorted, diag1_sorted, off2_sorted, diag2_sorted);
}

double dsk_distance_normalized(const ttk::DiagramType &D1,
                               const ttk::DiagramType &D2,
                               int L) {
  return std::sqrt(skot_distance_normalized(D1, D2, L));
}

void precompute_all_skot_sorted_vectors(
  const std::vector<ttk::DiagramType> &diagrams,
  int L,
  std::vector<std::vector<double>> &all_off_sorted,
  std::vector<std::vector<double>> &all_diag_sorted) {

  all_off_sorted.resize(diagrams.size());
  all_diag_sorted.resize(diagrams.size());

  for(std::size_t i = 0; i < diagrams.size(); ++i) {
    diagram_to_skot_sorted_vectors(
      diagrams[i], L, all_off_sorted[i], all_diag_sorted[i]);
  }
}

void precompute_all_skot_distances_from_empty(const std::vector<ttk::DiagramType> &diagrams,
                                              int L,
                                              std::vector<double> &distances, 
                                              std::vector<std::vector<double>> &all_off_sorted,
                                              std::vector<std::vector<double>> &all_diag_sorted) {
  distances.resize(diagrams.size());

  for(std::size_t i = 0; i < diagrams.size(); ++i) {
    distances[i] = dsk_distance_to_empty_from_precomputed_sorted(all_off_sorted[i], all_diag_sorted[i]);
  }
}

void precompute_all_w2_delta_sk_from_empty(const std::vector<ttk::DiagramType> &diagrams,
                                              int L,
                                              std::vector<double> &distances, 
                                              std::vector<std::vector<SkOffAtom>> &all_off_sorted,
                                              std::vector<std::vector<double>> &all_diag_sorted) {
  distances.resize(diagrams.size());

  for(std::size_t i = 0; i < diagrams.size(); ++i) {
    distances[i] = std::sqrt(w2_delta_sk_squared_to_empty_from_precomputed_sorted(all_off_sorted[i], all_diag_sorted[i]));
  }
}

void precompute_all_pot_distances_from_empty(
  const std::vector<std::vector<double>> &allV,
  int L,
  int choiceHilbertDistance,
  std::vector<double> &distances, int GLevel) {

  distances.resize(allV.size());

  if(choiceHilbertDistance == 1) {
    for(std::size_t i = 0; i < allV.size(); ++i) {
      double cost = 0.0;
      for(const double t : allV[i]) {
        cost += (t <= 0.5) ? t : (1.0 - t);
      }
      distances[i] = std::sqrt(cost);
    }
  } else if(choiceHilbertDistance == 3) {
    static const std::vector<double> empty;
    const int Pbits = std::min(60, L + 1);
    const int L4max = std::min(26, (Pbits + 1) / 2);

    for(std::size_t i = 0; i < allV.size(); ++i) {
      const int GnLevel = GLevel;
      const int L_eps   = L;
      const double cost = pot_figalli_1d_fast_L1_generatorsGn_cost_sorted(allV[i], empty, GnLevel, L_eps);
      distances[i] = std::sqrt(cost);
    }
  }
}

void diagram_to_v_tmid(const ttk::DiagramType &D,
                       int L,
                       std::vector<double> &out_v) {
  out_v.clear();
  out_v.reserve(D.size());

  auto clamp01_tt = [](long double v) -> long double {
    if(v < 0.0L) return 0.0L;
    if(v > 1.0L) return 1.0L;
    return v;
  };

  for(const auto &pt : D) {
    const long double b = (long double)pt.birth.sfValue;
    const long double d = (long double)pt.death.sfValue;
    if(!std::isfinite((double)b) || !std::isfinite((double)d)) continue;

    long double x = clamp01_tt(b);
    long double y = clamp01_tt(d);
    if(x > y) std::swap(x, y);
    if(y <= x) continue;

    out_v.push_back(hilbert_tmid(L, x, y));
  }

  std::sort(out_v.begin(), out_v.end());
}

std::vector<std::vector<double>>
precompute_all_v(const std::vector<ttk::DiagramType> &geodesic, int L) {
  std::vector<std::vector<double>> allV(geodesic.size());
  for(size_t i = 0; i < geodesic.size(); ++i) {
    diagram_to_v_tmid(geodesic[i], L, allV[i]);
  }
  return allV;
}

double pot_cost_on_SK_1d_from_sorted_tmid(const std::vector<double> &v1_sorted,
                                         const std::vector<double> &v2_sorted,
                                         int L,
                                         int choiceHilbertDistance, int GLevel) {
  const int Pbits = std::min(60, L + 1);
  const int L4max = std::min(26, (Pbits + 1) / 2);

  if(choiceHilbertDistance == 3) {
    const int GnLevel = GLevel;
    const int L_eps   = L;
    return pot_figalli_1d_fast_L1_generatorsGn_cost_sorted(v1_sorted, v2_sorted, GnLevel, L_eps);
  }
  
  return pot_figalli_1d_fast_L1_sorted_cost(v1_sorted, v2_sorted);
}


void diagram_to_w2_delta_sk_sorted_data(
  const ttk::DiagramType &D,
  int L,
  std::vector<SkOffAtom> &out_off_sorted,
  std::vector<double> &out_diag_sorted) {

  out_off_sorted.clear();
  out_diag_sorted.clear();

  out_off_sorted.reserve(D.size());
  out_diag_sorted.reserve(D.size());

  auto clamp01_tt = [](long double v) -> long double {
    if(v < 0.0L)
      return 0.0L;
    if(v > 1.0L)
      return 1.0L;
    return v;
  };

  for(const auto &pt : D) {
    const long double b = static_cast<long double>(pt.birth.sfValue);
    const long double d = static_cast<long double>(pt.death.sfValue);

    if(!std::isfinite(static_cast<double>(b)) ||
       !std::isfinite(static_cast<double>(d)))
      continue;

    long double x = clamp01_tt(b);
    long double y = clamp01_tt(d);
    if(x > y)
      std::swap(x, y);

    if(y <= x)
      continue;

    const double toff = sk_selector_value(L, x, y);
    out_off_sorted.push_back(
      SkOffAtom{toff, static_cast<double>(x), static_cast<double>(y)});

    const long double m = 0.5L * (x + y);
    out_diag_sorted.push_back(sk_selector_value(L, m, m));
  }

  std::sort(out_off_sorted.begin(), out_off_sorted.end(), sk_off_atom_less);
  std::sort(out_diag_sorted.begin(), out_diag_sorted.end());
}

double w2_delta_sk_squared_from_precomputed_sorted(
  const std::vector<SkOffAtom> &off1_sorted,
  const std::vector<double> &diag1_sorted,
  const std::vector<SkOffAtom> &off2_sorted,
  const std::vector<double> &diag2_sorted) {

  assert(off1_sorted.size() == diag1_sorted.size());
  assert(off2_sorted.size() == diag2_sorted.size());

  const std::size_t Nleft  = off1_sorted.size() + diag2_sorted.size();
  const std::size_t Nright = off2_sorted.size() + diag1_sorted.size();
  assert(Nleft == Nright);

  std::size_t i_off1 = 0, i_diag2 = 0;
  std::size_t i_off2 = 0, i_diag1 = 0;

  double cost_sq = 0.0;

  for(std::size_t k = 0; k < Nleft; ++k) {
    const SkAugAtom a = next_from_two_sorted_augmented(
      off1_sorted, i_off1, diag2_sorted, i_diag2);
    const SkAugAtom b = next_from_two_sorted_augmented(
      off2_sorted, i_off2, diag1_sorted, i_diag1);

    cost_sq += c_delta_sq(a, b);
  }

  return cost_sq;
}

double w2_delta_sk_from_precomputed_sorted(
  const std::vector<SkOffAtom> &off1_sorted,
  const std::vector<double> &diag1_sorted,
  const std::vector<SkOffAtom> &off2_sorted,
  const std::vector<double> &diag2_sorted) {

  return std::sqrt(
    w2_delta_sk_squared_from_precomputed_sorted(
      off1_sorted, diag1_sorted, off2_sorted, diag2_sorted));
}

double w2_delta_sk_squared_normalized(
  const ttk::DiagramType &D1,
  const ttk::DiagramType &D2,
  int L) {

  std::vector<SkOffAtom> off1_sorted, off2_sorted;
  std::vector<double> diag1_sorted, diag2_sorted;

  diagram_to_w2_delta_sk_sorted_data(D1, L, off1_sorted, diag1_sorted);
  diagram_to_w2_delta_sk_sorted_data(D2, L, off2_sorted, diag2_sorted);

  return w2_delta_sk_squared_from_precomputed_sorted(
    off1_sorted, diag1_sorted, off2_sorted, diag2_sorted);
}

double w2_delta_sk_normalized(
  const ttk::DiagramType &D1,
  const ttk::DiagramType &D2,
  int L) {

  return std::sqrt(w2_delta_sk_squared_normalized(D1, D2, L));
}

void precompute_all_w2_delta_sk_sorted_data(
  const std::vector<ttk::DiagramType> &diagrams,
  int L,
  std::vector<std::vector<SkOffAtom>> &all_off_sorted,
  std::vector<std::vector<double>> &all_diag_sorted) {

  all_off_sorted.resize(diagrams.size());
  all_diag_sorted.resize(diagrams.size());

  for(std::size_t i = 0; i < diagrams.size(); ++i) {
    diagram_to_w2_delta_sk_sorted_data(
      diagrams[i], L, all_off_sorted[i], all_diag_sorted[i]);
  }
}

double skot_distance_to_empty_normalized(
  const ttk::DiagramType &D,
  int L) {

  std::vector<double> off_sorted, diag_sorted;
  diagram_to_skot_sorted_vectors(D, L, off_sorted, diag_sorted);
  return skot_distance_to_empty_from_precomputed_sorted(off_sorted, diag_sorted);
}

double dsk_distance_to_empty_normalized(
  const ttk::DiagramType &D,
  int L) {

  return std::sqrt(skot_distance_to_empty_normalized(D, L));
}

double w2_delta_sk_squared_to_empty_normalized(
  const ttk::DiagramType &D,
  int L) {

  std::vector<SkOffAtom> off_sorted;
  std::vector<double> diag_sorted;
  diagram_to_w2_delta_sk_sorted_data(D, L, off_sorted, diag_sorted);

  return w2_delta_sk_squared_to_empty_from_precomputed_sorted(
    off_sorted, diag_sorted);
}

double w2_delta_sk_to_empty_normalized(
  const ttk::DiagramType &D,
  int L) {

  return std::sqrt(w2_delta_sk_squared_to_empty_normalized(D, L));
}

double skot_distance_to_empty_from_precomputed_sorted(
  const std::vector<double> &off_sorted,
  const std::vector<double> &diag_sorted) {

  assert(off_sorted.size() == diag_sorted.size());

  double cost = 0.0;
  for(std::size_t k = 0; k < off_sorted.size(); ++k) {
    cost += std::abs(off_sorted[k] - diag_sorted[k]);
  }
  return cost;
}

double dsk_distance_to_empty_from_precomputed_sorted(
  const std::vector<double> &off_sorted,
  const std::vector<double> &diag_sorted) {

  return std::sqrt(
    skot_distance_to_empty_from_precomputed_sorted(off_sorted, diag_sorted));
}

double w2_delta_sk_squared_to_empty_from_precomputed_sorted(
  const std::vector<SkOffAtom> &off_sorted,
  const std::vector<double> &diag_sorted) {

  assert(off_sorted.size() == diag_sorted.size());

  double cost_sq = 0.0;
  for(const auto &a : off_sorted) {
    const double pers = a.y - a.x;
    cost_sq += 0.5 * pers * pers;
  }
  return cost_sq;
}

double w2_delta_sk_to_empty_from_precomputed_sorted(
  const std::vector<SkOffAtom> &off_sorted,
  const std::vector<double> &diag_sorted) {

  return std::sqrt(
    w2_delta_sk_squared_to_empty_from_precomputed_sorted(
      off_sorted, diag_sorted));
}