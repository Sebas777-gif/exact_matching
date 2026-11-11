#pragma once

// Coarse arrangement A0 in the strip 0 <= s <= 1 of the (s, t) plane of
// slices y = s x + t: one dual segment per distinct input grade, plus the
// walls s = 0 and s = 1. Slopes s > 1 are covered by the same strip for the
// modules with x and y swapped (see swapped_module).

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <algorithm>
#include <cassert>
#include <functional>
#include <numeric>
#include <optional>
#include <set>
#include <utility>
#include <vector>

#include <phat/boundary_matrix.h>
#include <phat/persistence_pairs.h>
#include <phat/compute_persistence_pairs.h>

namespace matching {

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Segment_2 Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
typedef Traits_2::Point_2 Point_2;  // (s, t)
typedef Traits_2::X_monotone_curve_2 X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;

using Index = phat::index;
using IndexVec = std::vector<Index>;

struct Grade2 {
  FT x, y;
};

struct Bar {
  FT birth, death;
  bool infinite;
};

inline bool operator<(const Grade2& a, const Grade2& b) {
  if (a.x != b.x) return a.x < b.x;
  return a.y < b.y;
}

inline bool operator==(const Grade2& a, const Grade2& b) {
  return a.x == b.x && a.y == b.y;
}

// Finite presentation. columns is the boundary matrix; A0 uses only the grades.
struct Presentation {
  std::vector<Grade2> generators;
  std::vector<Grade2> relations;
  std::vector<std::vector<int> > columns;
};

// Reduce a list of generator indices mod 2: a generator listed twice
// cancels (g + g = 0), so equal indices are removed in pairs.
// The result is sorted.
inline void reduce_mod2(std::vector<int>& col) {
  std::sort(col.begin(), col.end());
  std::vector<int> out;
  out.reserve(col.size());
  for (int idx : col) {
    if (!out.empty() && out.back() == idx)
      out.pop_back();
    else
      out.push_back(idx);
  }
  col = std::move(out);
}

struct Module {
  std::vector<Presentation> summands;
};

inline Module as_module(Presentation P) {
  Module M;
  M.summands.push_back(std::move(P));
  return M;
}

// Dual of a grade v = (v_x, v_y): the slices through v, namely t = v_y - s v_x,
// clipped to the strip 0 <= s <= 1.
inline Segment_2 dual_segment_of_grade(const Grade2& v) {
  return Segment_2(Point_2(FT(0), v.y), Point_2(FT(1), v.y - v.x));
}

inline void collect_grades(const Module& M, const Module& N, std::set<Grade2>& grades) {
  grades.clear();
  auto add = [&](const Module& Mod) {
    for (const Presentation& P : Mod.summands) {
      for (const Grade2& g : P.generators) grades.insert(g);
      for (const Grade2& g : P.relations) grades.insert(g);
    }
  };
  add(M);
  add(N);
}

inline Module swapped_module(const Module& M) {
  Module swap_M;
  std::vector<Presentation> swap_summands;
  for (const Presentation& P : M.summands) {
    Presentation swap_P;
    std::vector<Grade2> gens;
    std::vector<Grade2> rels;
    for (const Grade2& g : P.generators) {
      Grade2 swap_g;
      swap_g.y = g.x;
      swap_g.x = g.y;
      gens.push_back(swap_g);
    }
    for (const Grade2& g : P.relations) {
      Grade2 swap_g;
      swap_g.y = g.x;
      swap_g.x = g.y;
      rels.push_back(swap_g);
    }
    swap_P.generators = gens;
    swap_P.relations = rels;
    swap_P.columns = P.columns;
    swap_summands.push_back(swap_P);
  }
  swap_M.summands = swap_summands;
  return swap_M;
}

// Inserts A0 restricted to the strip 0 <= s <= 1 and returns how many
// segments were inserted: one per distinct grade, plus the walls.
// Each wall s = 0 and s = 1 runs from the lowest to the highest endpoint of
// a dual segment on it. A wall whose endpoints all coincide is a single
// point, which is already a vertex, so it is left out.
inline std::size_t build_arrangement_A(Arrangement_2& arr, const Module& M, const Module& N) {
  std::set<Grade2> grades;
  collect_grades(M, N, grades);
  if (grades.empty()) return 0;

  std::vector<X_monotone_curve_2> curves;
  curves.reserve(grades.size() + 2);
  FT t0_min = grades.begin()->y, t0_max = t0_min;
  FT t1_min = grades.begin()->y - grades.begin()->x, t1_max = t1_min;
  for (const Grade2& g : grades) {
    Segment_2 seg = dual_segment_of_grade(g);
    curves.push_back(X_monotone_curve_2(seg));
    t0_min = std::min(t0_min, seg.source().y());
    t0_max = std::max(t0_max, seg.source().y());
    t1_min = std::min(t1_min, seg.target().y());
    t1_max = std::max(t1_max, seg.target().y());
  }
  if (t0_min < t0_max)
    curves.push_back(X_monotone_curve_2(Point_2(FT(0), t0_min), Point_2(FT(0), t0_max)));
  if (t1_min < t1_max)
    curves.push_back(X_monotone_curve_2(Point_2(FT(1), t1_min), Point_2(FT(1), t1_max)));
  CGAL::insert(arr, curves.begin(), curves.end());
  return curves.size();
}

inline FT weighted_push(const Grade2& p, const FT& s, const FT& t) {
  assert(s >= 0);
  if (p.y > s * p.x + t) { // strictly above the slice
    if (s >= 1)
      return (p.y - t) / s;
    return p.y;
  }
  if (s >= 1)
    return p.x;
  return s * p.x + t;
}

// adapted from Hera (persistence_module.hpp)

template<class T>
IndexVec get_sorted_indices(const std::vector<T>& values) { 
  IndexVec result(values.size());
  std::iota(result.begin(), result.end(), 0);
  std::sort(result.begin(), result.end(),
          [&values](Index a, Index b) { return values[a] < values[b]; });
  return result;
}

inline std::vector<FT> project_generators(const FT& s, const FT& t, const Presentation& P) {
  std::vector<FT> gen_values;
  gen_values.reserve(P.generators.size());
  for (const Grade2& g : P.generators) gen_values.push_back(weighted_push(g, s, t));
  return gen_values;
}

inline std::vector<FT> project_relations(const FT& s, const FT& t, const Presentation& P) {
  std::vector<FT> rel_values;
  rel_values.reserve(P.relations.size());
  for (const Grade2& r : P.relations) rel_values.push_back(weighted_push(r, s, t));
  return rel_values;
}

inline phat::boundary_matrix<> build_slice_matrix(const Presentation& P,
    const IndexVec& gen_order, const IndexVec& rel_order) {
  const Index n_gens = P.generators.size();
  const Index n_rels = P.relations.size();
  phat::boundary_matrix<> phat_matrix;

  IndexVec gen_rank(n_gens);
  for (Index k = 0; k < n_gens; k++) gen_rank[gen_order[k]] = k;

  phat_matrix.set_num_cols(n_gens + n_rels);
  for (Index k = 0; k < n_gens; k++) phat_matrix.set_dim(k, 0);
  for (Index i = 0; i < n_rels; i++) {
    const auto& col = P.columns[rel_order[i]];
    IndexVec current_rel(col.begin(), col.end());
    for (Index& j : current_rel) j = gen_rank[j];
    std::sort(current_rel.begin(), current_rel.end());
    phat_matrix.set_dim(n_gens + i, 1);
    phat_matrix.set_col(n_gens + i, current_rel);
  }
  return phat_matrix;
}

inline std::vector<Bar> slice_barcode(const Presentation& P, const FT& s, const FT& t) {
  std::vector<FT> proj_gen = project_generators(s, t, P);
  std::vector<FT> proj_rel = project_relations(s, t, P);
  IndexVec gen_order = get_sorted_indices(proj_gen);
  IndexVec rel_order = get_sorted_indices(proj_rel);
  phat::boundary_matrix<> phat_matrix = build_slice_matrix(P, gen_order, rel_order);
  phat::persistence_pairs pairs;
  phat::compute_persistence_pairs<phat::twist_reduction>(pairs, phat_matrix);

  Index n_gens = P.generators.size();
  std::vector<Bar> bars;
  // each unmarked i < n_gens will be an infinite bar starting at gen_values[gen_order[i]]
  std::vector<bool> is_finite(n_gens, false); 
  for (Index i = 0; i < pairs.get_num_pairs(); i++) {
    auto p = pairs.get_pair(i);
    Index k = p.first;
    Index j = p.second;
    is_finite[k] = true;
    Bar b;
    b.birth = proj_gen[gen_order[k]];
    b.death = proj_rel[rel_order[j - n_gens]];
    b.infinite = false;
    if (b.birth != b.death) bars.push_back(b);
  }
  for (Index i = 0; i < n_gens; i++)
    if (!is_finite[i]) {
      Bar b;
      b.birth = proj_gen[gen_order[i]];
      b.infinite = true;
      bars.push_back(b);
  }
  return bars;
}

namespace detail {

// Is there a matching of cost at most c between the finite barcodes A and B,
// where a bar may also be matched to the diagonal? Costs are given as ranks
// into the sorted list of candidate values, so c is a rank too:
//   pair_rank[i][j]  rank of the cost of matching A[i] with B[j]
//   diag_rank_A[i]   rank of the cost of matching A[i] to the diagonal
//   diag_rank_B[j]   rank of the cost of matching B[j] to the diagonal
//
// This is a perfect matching problem on a bipartite graph with n + m
// vertices per side:
//   left:  A[0], ..., A[n-1], then a diagonal copy of each B[j] at n + j
//   right: B[0], ..., B[m-1], then a diagonal copy of each A[i] at m + i
// A[i] can go to B[j] or to its own diagonal copy, B[j] can go to A[i] or to
// its own diagonal copy, and diagonal copies can always be matched to each
// other (cost 0). It is solved with augmenting paths (Kuhn's algorithm).
inline bool bottleneck_feasible(std::size_t c,
    const std::vector<std::vector<std::size_t> >& pair_rank,
    const std::vector<std::size_t>& diag_rank_A,
    const std::vector<std::size_t>& diag_rank_B) {
  const std::size_t n = diag_rank_A.size();
  const std::size_t m = diag_rank_B.size();
  const std::size_t N = n + m;

  auto edge = [&](std::size_t l, std::size_t r) {
    if (l < n && r < m) return pair_rank[l][r] <= c;   // A[l] with B[r]
    if (l < n) return r - m == l && diag_rank_A[l] <= c; // A[l] to diagonal
    if (r < m) return l - n == r && diag_rank_B[r] <= c; // B[r] to diagonal
    return true;                                         // diagonal with diagonal
  };

  const std::size_t none = N;
  std::vector<std::size_t> match_of_right(N, none);
  std::vector<char> visited(N);

  std::function<bool(std::size_t)> augment = [&](std::size_t l) {
    for (std::size_t r = 0; r < N; r++) {
      if (visited[r] || !edge(l, r)) continue;
      visited[r] = 1;
      if (match_of_right[r] == none || augment(match_of_right[r])) {
        match_of_right[r] = l;
        return true;
      }
    }
    return false;
  };

  for (std::size_t l = 0; l < N; l++) {
    std::fill(visited.begin(), visited.end(), 0);
    if (!augment(l)) return false;
  }
  return true;
}

}  // namespace detail

// Exact bottleneck distance between two barcodes (L-infinity cost, a bar may
// be matched to the diagonal at cost (death - birth) / 2).
// Returns std::nullopt, meaning infinity, if the numbers of infinite bars differ.
//
// Infinite bars are matched to each other in order of birth. For the finite
// bars the distance is one of the candidate values (a pair cost or a diagonal
// cost); a binary search over the sorted candidates finds the smallest one
// for which a matching exists.
inline std::optional<FT> bottleneck_dist(const std::vector<Bar>& b1, const std::vector<Bar>& b2) {
  std::vector<Bar> A, B;
  std::vector<FT> inf_A, inf_B;
  for (const Bar& b : b1) {
    if (b.infinite) inf_A.push_back(b.birth);
    else A.push_back(b);
  }
  for (const Bar& b : b2) {
    if (b.infinite) inf_B.push_back(b.birth);
    else B.push_back(b);
  }

  if (inf_A.size() != inf_B.size()) return std::nullopt;
  std::sort(inf_A.begin(), inf_A.end());
  std::sort(inf_B.begin(), inf_B.end());
  FT result(0);
  for (std::size_t k = 0; k < inf_A.size(); k++)
    result = std::max(result, CGAL::abs(inf_A[k] - inf_B[k]));

  const std::size_t n = A.size();
  const std::size_t m = B.size();
  if (n + m == 0) return result;

  auto diag_cost = [](const Bar& b) { return (b.death - b.birth) / 2; };
  auto pair_cost = [](const Bar& a, const Bar& b) {
    return std::max(CGAL::abs(a.birth - b.birth), CGAL::abs(a.death - b.death));
  };

  std::vector<FT> candidates;
  candidates.reserve(n * m + n + m);
  for (const Bar& a : A) candidates.push_back(diag_cost(a));
  for (const Bar& b : B) candidates.push_back(diag_cost(b));
  for (const Bar& a : A)
    for (const Bar& b : B) candidates.push_back(pair_cost(a, b));
  std::sort(candidates.begin(), candidates.end());
  candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());

  // Replace every cost by its rank, so the matching only compares integers.
  auto rank = [&](const FT& x) -> std::size_t {
    return std::lower_bound(candidates.begin(), candidates.end(), x) - candidates.begin();
  };
  std::vector<std::size_t> diag_rank_A(n), diag_rank_B(m);
  std::vector<std::vector<std::size_t> > pair_rank(n, std::vector<std::size_t>(m));
  for (std::size_t i = 0; i < n; i++) diag_rank_A[i] = rank(diag_cost(A[i]));
  for (std::size_t j = 0; j < m; j++) diag_rank_B[j] = rank(diag_cost(B[j]));
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < m; j++) pair_rank[i][j] = rank(pair_cost(A[i], B[j]));

  // The largest candidate always works: match every bar to the diagonal.
  std::size_t lo = 0, hi = candidates.size() - 1;
  while (lo < hi) {
    std::size_t mid = (lo + hi) / 2;
    if (detail::bottleneck_feasible(mid, pair_rank, diag_rank_A, diag_rank_B)) hi = mid;
    else lo = mid + 1;
  }
  return std::max(result, candidates[lo]);
}

inline std::vector<Bar> module_barcode(const Module& M, const FT& s, const FT& t) {
  std::vector<Bar> bars;
  for (const Presentation& P : M.summands) {
    std::vector<Bar> b = slice_barcode(P, s, t);
    bars.insert(bars.end(), b.begin(), b.end());
  }
  return bars;
}

inline std::optional<FT> F(const Module& M, const Module& N, const FT& s, const FT& t) {
  return bottleneck_dist(module_barcode(M, s, t), module_barcode(N, s, t));
}

}  // namespace matching
