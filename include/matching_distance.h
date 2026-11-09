#pragma once

// Coarse arrangement A0 in the (s, t) plane of slices y = s x + t.
// One dual ray per distinct input grade, plus the walls s = 0 and s = 1.

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <algorithm>
#include <cassert>
#include <numeric>
#include <set>
#include <utility>
#include <vector>

#include <phat/boundary_matrix.h>
#include <phat/persistence_pairs.h>
#include <phat/compute_persistence_pairs.h>

namespace matching {

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef CGAL::Arr_linear_traits_2<Kernel> Traits_2;
typedef Traits_2::Point_2 Point_2;  // (s, t)
typedef Traits_2::Line_2 Line_2;    // s = 0 and s = 1
typedef Traits_2::X_monotone_curve_2 X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;

using Index = phat::index;
using IndexVec = std::vector<Index>;

struct Grade2 {
  FT x, y;
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

// Dual of a grade v = (v_x, v_y): the slices through v, namely t = v_y - s v_x
inline Line_2 dual_line_of_grade(const Grade2& v) {
  return Line_2(v.x, FT(1), -v.y);
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

// Inserts A0 and returns how many lines were inserted.
inline std::size_t build_arrangement_A(Arrangement_2& arr, const Module& M, const Module& N) {
  std::set<Grade2> grades;
  collect_grades(M, N, grades);

  std::vector<X_monotone_curve_2> curves;
  curves.reserve(grades.size() + 2);
  for (const Grade2& g : grades)
    curves.push_back(X_monotone_curve_2(dual_line_of_grade(g)));
  curves.push_back(X_monotone_curve_2(Line_2(Point_2(FT(1), FT(0)), Point_2(FT(1), FT(1)))));
  curves.push_back(X_monotone_curve_2(Line_2(Point_2(FT(0), FT(0)), Point_2(FT(0), FT(1)))));
  CGAL::insert(arr, curves.begin(), curves.end());
  return curves.size();
}

inline FT weighted_push(const Grade2& p, const FT& s, const FT& t) {
  assert(s > 0);
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

}  // namespace matching
