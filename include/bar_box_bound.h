#pragma once

// The bar-box bound: an upper bound on F over a convex polygon of slices,
// without evaluating F. See notes/bar_box_bound/bar_box_bound.pdf.
//
// For a polygon with corners c_i, every grade's push takes values in the
// range [min_i push(c_i), max_i push(c_i)] on the whole polygon. Within each
// summand, overlapping ranges of generators merge into generator classes, and
// likewise for relations. The number of bars with birth in a given generator
// class and death in a given relation class is the same at every slice of the
// polygon (pivot lemma). So the barcode of a module on the polygon is
// described by a fixed list of boxes J x K, and a bottleneck matching between
// the two lists of boxes bounds F on the whole polygon.
//
// The polygon must lie in the strip 0 <= s <= 1 (flat pushes; the steep
// strip uses the swapped modules).
//
// Arithmetic: every computation runs first in interval arithmetic. Whenever
// interval arithmetic cannot decide a comparison, the conservative answer is
// taken (classes merge, a box is kept, an edge of the matching is dropped, a
// box counts as forced), which only makes the bound weaker, never wrong. The
// green test then repeats the computation exactly if that happened, so ties
// such as a bound exactly equal to Delta are decided correctly.

#include <CGAL/Interval_nt.h>

#include <algorithm>
#include <functional>
#include <limits>
#include <numeric>
#include <set>
#include <tuple>
#include <vector>

#include <matching_distance.h>

namespace matching {

typedef CGAL::Interval_nt<true> Interval;

// One summand, as the bar-box bound needs it: its grades and its matrix
// (columns[r] = generator indices of relation r, reduced mod 2).
struct BarBoxSummand {
  std::vector<Grade2> gens, rels;
  std::vector<std::vector<int> > columns;
};
typedef std::vector<BarBoxSummand> BarBoxModule;

// Summands without generators have no bars and are left out.
inline BarBoxModule bar_box_module(const Module& M) {
  BarBoxModule out;
  for (const Presentation& P : M.summands) {
    if (P.generators.empty()) continue;
    out.push_back(BarBoxSummand{P.generators, P.relations, P.columns});
  }
  return out;
}

// --- Arithmetic policies ------------------------------------------------------
//
// le(a, b) and lt(a, b) return true only if the comparison certainly holds.
// IntervalArithmetic records in `uncertain` that it could not decide.

struct IntervalArithmetic {
  typedef Interval T;
  bool uncertain = false;

  static T from(const FT& x) { return T(CGAL::to_interval(x)); }
  static T zero() { return T(0); }
  static T half(const T& x) { return x * T(0.5); }
  static T max(const T& a, const T& b) {
    return T(std::max(a.inf(), b.inf()), std::max(a.sup(), b.sup()));
  }
  static T min(const T& a, const T& b) {
    return T(std::min(a.inf(), b.inf()), std::min(a.sup(), b.sup()));
  }
  static bool sort_before(const T& a, const T& b) { return a.inf() < b.inf(); }
  // push(v) = max(v_y, s v_x + t)
  static T push(const Grade2& v, const T& s, const T& t) { return max(from(v.y), s * from(v.x) + t); }

  bool le(const T& a, const T& b) {
    if (a.sup() <= b.inf()) return true;
    if (a.inf() > b.sup()) return false;
    uncertain = true;
    return false;
  }
  bool lt(const T& a, const T& b) {
    if (a.sup() < b.inf()) return true;
    if (a.inf() >= b.sup()) return false;
    uncertain = true;
    return false;
  }
};

struct ExactArithmetic {
  typedef FT T;
  bool uncertain = false;  // never set

  static T from(const FT& x) { return x; }
  static T zero() { return FT(0); }
  static T half(const T& x) { return x / 2; }
  static T max(const T& a, const T& b) { return std::max(a, b); }
  static T min(const T& a, const T& b) { return std::min(a, b); }
  static bool sort_before(const T& a, const T& b) { return a < b; }
  static T push(const Grade2& v, const T& s, const T& t) { return weighted_push(v, s, t); }

  bool le(const T& a, const T& b) { return a <= b; }
  bool lt(const T& a, const T& b) { return a < b; }
};

// --- Boxes ------------------------------------------------------------------------

// A box J x K of a pivot pair (birth in J, death in K), or J x {infinity}
// for an unpaired generator (essential bar).
template <class T>
struct BarBox {
  T jlo, jhi, klo, khi;
  bool essential;
  std::size_t summand;
  int gen_class, rel_class;  // rel_class is -1 for essential boxes
};

// The boxes of one module on a polygon, with multiplicity, and the members of
// every class (for the local phase).
template <class T>
struct BarBoxes {
  std::vector<BarBox<T> > boxes;
  // gen_members[q][a]: generator indices of class a of summand q; likewise rel_members.
  std::vector<std::vector<std::vector<int> > > gen_members, rel_members;
};

namespace detail {

// Groups grades into classes: the connected components of the union of their
// ranges [lo, hi]. Returns the class of every grade and the hull of every class.
template <class A>
void bar_box_classes(const std::vector<typename A::T>& lo, const std::vector<typename A::T>& hi,
                     A& arith, std::vector<int>& class_of,
                     std::vector<std::pair<typename A::T, typename A::T> >& hull) {
  const std::size_t n = lo.size();
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](int a, int b) { return A::sort_before(lo[a], lo[b]); });
  class_of.assign(n, -1);
  hull.clear();
  for (int g : order) {
    // A new class starts only if it is certainly separated from the last one.
    if (hull.empty() || arith.lt(hull.back().second, lo[g])) {
      hull.push_back({lo[g], hi[g]});
    } else {
      hull.back().first = A::min(hull.back().first, lo[g]);
      hull.back().second = A::max(hull.back().second, hi[g]);
    }
    class_of[g] = static_cast<int>(hull.size()) - 1;
  }
}

// The pivot pairs of the column reduction with generators and relations in the
// given orders: (generator, relation) for every pair and (generator, -1) for
// every unpaired generator. The pivot of a column is its latest generator.
inline void pivot_pairs(const BarBoxSummand& q, const std::vector<int>& gen_order,
                        const std::vector<int>& rel_order, std::vector<std::pair<int, int> >& out) {
  const std::size_t n = q.gens.size();
  std::vector<int> rank(n);
  for (std::size_t i = 0; i < n; i++) rank[gen_order[i]] = static_cast<int>(i);
  std::vector<std::vector<int> > reduced;
  std::vector<int> owner(n, -1);  // owner[p]: reduced column with pivot p
  std::vector<char> paired(n, 0);
  std::vector<int> col, tmp;
  out.clear();
  for (int r : rel_order) {
    col.clear();
    for (int g : q.columns[r]) col.push_back(rank[g]);
    std::sort(col.begin(), col.end());
    while (!col.empty() && owner[col.back()] >= 0) {
      const std::vector<int>& other = reduced[owner[col.back()]];
      tmp.clear();
      std::set_symmetric_difference(col.begin(), col.end(), other.begin(), other.end(),
                                    std::back_inserter(tmp));
      col.swap(tmp);
    }
    if (col.empty()) continue;
    const int p = col.back();
    owner[p] = static_cast<int>(reduced.size());
    reduced.push_back(col);
    paired[p] = 1;
    out.push_back({gen_order[p], r});
  }
  for (std::size_t i = 0; i < n; i++)
    if (!paired[i]) out.push_back({gen_order[i], -1});
}

}  // namespace detail

// The boxes of module `mod` on the polygon with corners (s[i], t[i]).
template <class A>
BarBoxes<typename A::T> bar_boxes(const BarBoxModule& mod, const std::vector<typename A::T>& s,
                                  const std::vector<typename A::T>& t, A& arith) {
  typedef typename A::T T;
  BarBoxes<T> out;
  std::vector<T> glo, ghi, rlo, rhi;
  std::vector<int> gclass, rclass;
  std::vector<std::pair<T, T> > ghull, rhull;
  std::vector<std::pair<int, int> > pairs;

  auto ranges = [&](const std::vector<Grade2>& grades, std::vector<T>& lo, std::vector<T>& hi) {
    lo.clear();
    hi.clear();
    for (const Grade2& v : grades) {
      T l = A::push(v, s[0], t[0]), h = l;
      for (std::size_t i = 1; i < s.size(); i++) {
        const T p = A::push(v, s[i], t[i]);
        l = A::min(l, p);
        h = A::max(h, p);
      }
      lo.push_back(l);
      hi.push_back(h);
    }
  };

  for (std::size_t qi = 0; qi < mod.size(); qi++) {
    const BarBoxSummand& q = mod[qi];
    ranges(q.gens, glo, ghi);
    ranges(q.rels, rlo, rhi);
    detail::bar_box_classes(glo, ghi, arith, gclass, ghull);
    detail::bar_box_classes(rlo, rhi, arith, rclass, rhull);

    // Any order that sorts by class gives the same counts per pair of classes.
    std::vector<int> gen_order(q.gens.size()), rel_order(q.rels.size());
    std::iota(gen_order.begin(), gen_order.end(), 0);
    std::iota(rel_order.begin(), rel_order.end(), 0);
    std::stable_sort(gen_order.begin(), gen_order.end(),
                     [&](int a, int b) { return gclass[a] < gclass[b]; });
    std::stable_sort(rel_order.begin(), rel_order.end(),
                     [&](int a, int b) { return rclass[a] < rclass[b]; });
    detail::pivot_pairs(q, gen_order, rel_order, pairs);

    for (const auto& p : pairs) {
      const auto& J = ghull[gclass[p.first]];
      if (p.second < 0) {
        out.boxes.push_back(BarBox<T>{J.first, J.second, T(), T(), true, qi, gclass[p.first], -1});
        continue;
      }
      const auto& K = rhull[rclass[p.second]];
      // A box with max K <= min J only holds bars of length zero.
      if (arith.le(K.second, J.first)) continue;
      out.boxes.push_back(BarBox<T>{J.first, J.second, K.first, K.second, false, qi,
                                    gclass[p.first], rclass[p.second]});
    }

    out.gen_members.emplace_back(ghull.size());
    out.rel_members.emplace_back(rhull.size());
    for (std::size_t g = 0; g < q.gens.size(); g++) out.gen_members.back()[gclass[g]].push_back(g);
    for (std::size_t r = 0; r < q.rels.size(); r++) out.rel_members.back()[rclass[r]].push_back(r);
  }
  return out;
}

// Largest cost of sending a bar of box x to the diagonal: (max K - min J) / 2.
template <class A>
typename A::T bar_box_diag(const BarBox<typename A::T>& x) {
  return A::max(A::zero(), A::half(x.khi - x.jlo));
}

// Largest L-infinity distance between a point of box x and a point of box y.
// Only defined if both are finite or both are essential.
template <class A>
typename A::T bar_box_cost(const BarBox<typename A::T>& x, const BarBox<typename A::T>& y) {
  typename A::T c = A::max(A::zero(), A::max(x.jhi - y.jlo, y.jhi - x.jlo));
  if (!x.essential) c = A::max(c, A::max(x.khi - y.klo, y.khi - x.klo));
  return c;
}

// Whether the boxes of X and Y can be matched with every matched pair costing
// at most `threshold` and every unmatched box having diag <= threshold.
//
// A box is forced if it is essential or its diag exceeds the threshold. By
// the Mendelsohn-Dulmage theorem, such a matching exists if and only if one
// matching covers the forced boxes of X and one covers those of Y. Both are
// found with augmenting paths that start only at forced boxes.
template <class A>
bool bar_box_feasible(const std::vector<BarBox<typename A::T> >& X,
                      const std::vector<BarBox<typename A::T> >& Y,
                      const typename A::T& threshold, A& arith) {
  const std::size_t nx = X.size(), ny = Y.size();
  std::vector<signed char> edge(nx * ny, -1);  // cached allowed(x, y)
  auto allowed = [&](std::size_t x, std::size_t y) {
    signed char& e = edge[x * ny + y];
    if (e < 0)
      e = X[x].essential == Y[y].essential && arith.le(bar_box_cost<A>(X[x], Y[y]), threshold);
    return e == 1;
  };
  auto forced = [&](const BarBox<typename A::T>& b) {
    return b.essential || !arith.le(bar_box_diag<A>(b), threshold);
  };

  // Matches every forced box of one side into the other side.
  auto cover = [&](bool from_x) {
    const std::size_t n_from = from_x ? nx : ny, n_to = from_x ? ny : nx;
    std::vector<long> partner(n_to, -1);
    std::vector<char> seen(n_to);
    std::function<bool(std::size_t)> augment = [&](std::size_t u) {
      for (std::size_t v = 0; v < n_to; v++) {
        if (seen[v] || !(from_x ? allowed(u, v) : allowed(v, u))) continue;
        seen[v] = 1;
        if (partner[v] < 0 || augment(static_cast<std::size_t>(partner[v]))) {
          partner[v] = static_cast<long>(u);
          return true;
        }
      }
      return false;
    };
    for (std::size_t u = 0; u < n_from; u++) {
      if (!forced(from_x ? X[u] : Y[u])) continue;
      std::fill(seen.begin(), seen.end(), 0);
      if (!augment(u)) return false;
    }
    return true;
  };
  return cover(true) && cover(false);
}

// Corners of a polygon in the arithmetic of A.
template <class A>
void bar_box_corners(const std::vector<Point_2>& corners, std::vector<typename A::T>& s,
                     std::vector<typename A::T>& t) {
  s.clear();
  t.clear();
  for (const Point_2& c : corners) {
    s.push_back(A::from(c.x()));
    t.push_back(A::from(c.y()));
  }
}

// Whether the bar-box bound shows F <= delta on the whole polygon. Decided
// in interval arithmetic, and exactly if interval arithmetic was unsure.
inline bool bar_box_green(const BarBoxModule& M, const BarBoxModule& N,
                          const std::vector<Point_2>& corners, const FT& delta) {
  {
    IntervalArithmetic arith;
    std::vector<Interval> s, t;
    bar_box_corners<IntervalArithmetic>(corners, s, t);
    const auto X = bar_boxes(M, s, t, arith), Y = bar_boxes(N, s, t, arith);
    if (bar_box_feasible(X.boxes, Y.boxes, IntervalArithmetic::from(delta), arith)) return true;
    if (!arith.uncertain) return false;
  }
  ExactArithmetic arith;
  std::vector<FT> s, t;
  bar_box_corners<ExactArithmetic>(corners, s, t);
  const auto X = bar_boxes(M, s, t, arith), Y = bar_boxes(N, s, t, arith);
  return bar_box_feasible(X.boxes, Y.boxes, delta, arith);
}

// An upper bound on F over the polygon: the smallest candidate value at which
// the interval boxes can be matched. Infinite if no matching exists.
inline double bar_box_upper(const std::vector<BarBox<Interval> >& X,
                            const std::vector<BarBox<Interval> >& Y) {
  std::vector<double> values{0};
  for (const auto& x : X)
    if (!x.essential) values.push_back(bar_box_diag<IntervalArithmetic>(x).sup());
  for (const auto& y : Y)
    if (!y.essential) values.push_back(bar_box_diag<IntervalArithmetic>(y).sup());
  for (const auto& x : X)
    for (const auto& y : Y)
      if (x.essential == y.essential) values.push_back(bar_box_cost<IntervalArithmetic>(x, y).sup());
  std::sort(values.begin(), values.end());
  values.erase(std::unique(values.begin(), values.end()), values.end());
  IntervalArithmetic arith;
  auto feasible = [&](double v) { return bar_box_feasible(X, Y, Interval(v), arith); };
  if (!feasible(values.back())) return std::numeric_limits<double>::infinity();
  std::size_t lo = 0, hi = values.size() - 1;
  while (lo < hi) {
    const std::size_t mid = (lo + hi) / 2;
    if (feasible(values[mid])) hi = mid;
    else lo = mid + 1;
  }
  return values[lo];
}

}  // namespace matching
