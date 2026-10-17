#pragma once

#include <algorithm>
#include <functional>
#include <iostream>
#include <set>
#include <utility>
#include <vector>

#include <phat/algorithms/standard_reduction.h>
#include <phat/boundary_matrix.h>
#include <phat/compute_persistence_pairs.h>
#include <phat/persistence_pairs.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/intersections.h>
#include <CGAL/number_utils.h>

namespace matching {

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef CGAL::Arr_linear_traits_2<Kernel> Traits_2;
typedef Traits_2::Point_2 Point_2;
typedef Traits_2::Line_2 Line_2;
typedef Traits_2::Ray_2 Ray_2;
typedef Traits_2::Vector_2 Vector_2;
typedef Traits_2::X_monotone_curve_2 X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;

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

struct Presentation {
  std::vector<Grade2> generators;
  std::vector<Grade2> relations;
  std::vector<std::vector<int> > columns;
};

inline FT abs_ft(const FT& v) { return v < 0 ? -v : v; }

inline FT max_ft(const FT& a, const FT& b) { return a < b ? b : a; }

// Exact FT has no +∞. Negative bottleneck / F means an infinite value
// (slice diagrams with different numbers of infinite bars).
inline FT infinite_F() { return FT(-1); }
inline bool is_infinite_F(const FT& d) { return d < FT(0); }

enum DeltaType { TYPE_I, TYPE_II, TYPE_IIIA, TYPE_IIIB };

inline bool point_on_or_above_slice(const Grade2& p, const FT& s, const FT& t) {
  return p.y >= s * p.x + t;
}

inline DeltaType delta_type(const Grade2& p, const Grade2& q, const FT& s, const FT& t) {
  const bool p_above = point_on_or_above_slice(p, s, t);
  const bool q_above = point_on_or_above_slice(q, s, t);
  if (p_above && q_above) return TYPE_I;
  if (!p_above && !q_above) return TYPE_II;
  if (p_above) return TYPE_IIIA;
  return TYPE_IIIB;
}

// Weighted 1-parameter grade of a 2-parameter grade, equal to
// w(l) · lambda^-1(push(p, l)).  Bottleneck of these grades is F itself.
inline FT weighted_grade(const Grade2& p, const FT& s, const FT& t) {
  const bool above = point_on_or_above_slice(p, s, t);
  if (s >= FT(1)) {
    return above ? (p.y - t) / s : p.x;
  }
  return above ? (p.y - t) : s * p.x;
}

struct Bar {
  FT birth;
  FT death;
  bool infinite;
  Grade2 birth_pt;
  Grade2 death_pt;
  bool has_death_pt;
};

inline std::vector<Bar> slice_diagram(const Presentation& P, const FT& s, const FT& t) {
  const int n_gen = static_cast<int>(P.generators.size());
  const int n_rel = static_cast<int>(P.relations.size());
  std::vector<Bar> bars;
  if (n_gen == 0) return bars;

  std::vector<FT> gen_val(n_gen), rel_val(n_rel);
  for (int i = 0; i < n_gen; ++i) gen_val[i] = weighted_grade(P.generators[i], s, t);
  for (int i = 0; i < n_rel; ++i) rel_val[i] = weighted_grade(P.relations[i], s, t);

  std::vector<int> gen_order(n_gen), rel_order(n_rel);
  for (int i = 0; i < n_gen; ++i) gen_order[i] = i;
  for (int i = 0; i < n_rel; ++i) rel_order[i] = i;
  std::sort(gen_order.begin(), gen_order.end(),
            [&](int a, int b) {
              if (gen_val[a] != gen_val[b]) return gen_val[a] < gen_val[b];
              return a < b;
            });
  std::sort(rel_order.begin(), rel_order.end(),
            [&](int a, int b) {
              if (rel_val[a] != rel_val[b]) return rel_val[a] < rel_val[b];
              return a < b;
            });

  std::vector<int> gen_rank(n_gen);
  for (int r = 0; r < n_gen; ++r) gen_rank[gen_order[r]] = r;

  const int n_cols = n_gen + n_rel;
  phat::boundary_matrix<> matrix;
  matrix.set_num_cols(n_cols);
  for (int i = 0; i < n_gen; ++i) {
    matrix.set_dim(i, 0);
    matrix.set_col(i, phat::column());
  }
  for (int j = 0; j < n_rel; ++j) {
    const int col = n_gen + j;
    matrix.set_dim(col, 1);
    std::vector<phat::index> faces;
    const int orig = rel_order[j];
    for (int g : P.columns[orig]) {
      if (g < 0 || g >= n_gen) continue;
      faces.push_back(static_cast<phat::index>(gen_rank[g]));
    }
    std::sort(faces.begin(), faces.end());
    faces.erase(std::unique(faces.begin(), faces.end()), faces.end());
    matrix.set_col(col, faces);
  }

  phat::persistence_pairs pairs;
  phat::compute_persistence_pairs<phat::standard_reduction>(pairs, matrix);

  std::vector<char> gen_paired(n_gen, 0);
  for (phat::index k = 0; k < pairs.get_num_pairs(); ++k) {
    const std::pair<phat::index, phat::index> pr = pairs.get_pair(k);
    const int birth_col = static_cast<int>(pr.first);
    const int death_col = static_cast<int>(pr.second);
    if (birth_col < 0 || birth_col >= n_gen) continue;
    if (death_col < n_gen) continue;
    gen_paired[birth_col] = 1;
    const int gen_id = gen_order[birth_col];
    const int rel_sorted = death_col - n_gen;
    const int rel_id = rel_order[rel_sorted];
    Bar b;
    b.birth = gen_val[gen_id];
    b.death = rel_val[rel_id];
    b.infinite = false;
    b.birth_pt = P.generators[gen_id];
    b.death_pt = P.relations[rel_id];
    b.has_death_pt = true;
    if (b.birth < b.death) bars.push_back(b);
  }
  for (int r = 0; r < n_gen; ++r) {
    if (gen_paired[r]) continue;
    const int gen_id = gen_order[r];
    Bar b;
    b.birth = gen_val[gen_id];
    b.death = FT(0);
    b.infinite = true;
    b.birth_pt = P.generators[gen_id];
    b.has_death_pt = false;
    bars.push_back(b);
  }
  return bars;
}

struct Witness {
  Grade2 p, q;
  bool same_module;
  DeltaType type;
  bool valid;
};

inline FT linf_bars(const Bar& a, const Bar& b) {
  if (a.infinite && b.infinite) return abs_ft(a.birth - b.birth);
  if (a.infinite || b.infinite) return FT(-1); // infeasible
  return max_ft(abs_ft(a.birth - b.birth), abs_ft(a.death - b.death));
}

inline FT diag_cost(const Bar& a) {
  if (a.infinite) return FT(-1);
  return (a.death - a.birth) / FT(2);
}

inline bool can_match(const Bar& a, const Bar& b, const FT& eps) {
  const FT c = linf_bars(a, b);
  if (c < FT(0)) return false;
  return c <= eps;
}

inline bool can_diag(const Bar& a, const FT& eps) {
  const FT c = diag_cost(a);
  if (c < FT(0)) return false;
  return c <= eps;
}

inline bool matching_feasible(const std::vector<Bar>& A,
                              const std::vector<Bar>& B,
                              const FT& eps) {
  const int nA = static_cast<int>(A.size());
  const int nB = static_cast<int>(B.size());
  std::vector<char> mustA(nA), mustB(nB);
  std::vector<std::vector<int> > adjA(nA), adjB(nB);
  for (int i = 0; i < nA; ++i) {
    mustA[i] = !can_diag(A[i], eps);
    for (int j = 0; j < nB; ++j) {
      if (can_match(A[i], B[j], eps)) {
        adjA[i].push_back(j);
        adjB[j].push_back(i);
      }
    }
  }
  for (int j = 0; j < nB; ++j) mustB[j] = !can_diag(B[j], eps);

  std::vector<int> ma(nA, -1), mb(nB, -1);
  std::vector<char> seenA, seenB;

  std::function<bool(int)> dfs_left = [&](int u) -> bool {
    for (int v : adjA[u]) {
      if (seenB[v]) continue;
      seenB[v] = 1;
      if (mb[v] == -1 || dfs_left(mb[v])) {
        ma[u] = v;
        mb[v] = u;
        return true;
      }
    }
    return false;
  };

  std::function<bool(int)> dfs_right = [&](int v) -> bool {
    for (int u : adjB[v]) {
      if (seenA[u]) continue;
      seenA[u] = 1;
      if (ma[u] == -1 || dfs_right(ma[u])) {
        ma[u] = v;
        mb[v] = u;
        return true;
      }
    }
    return false;
  };

  for (int u = 0; u < nA; ++u) {
    if (!mustA[u]) continue;
    seenB.assign(nB, 0);
    if (!dfs_left(u)) return false;
  }
  for (int v = 0; v < nB; ++v) {
    if (!mustB[v] || mb[v] != -1) continue;
    seenA.assign(nA, 0);
    if (!dfs_right(v)) return false;
  }
  return true;
}

inline FT bottleneck(const std::vector<Bar>& A, const std::vector<Bar>& B, Witness* witness) {
  int infA = 0, infB = 0;
  for (const Bar& b : A) if (b.infinite) ++infA;
  for (const Bar& b : B) if (b.infinite) ++infB;
  if (infA != infB) {
    if (witness) witness->valid = false;
    return infinite_F();
  }

  std::vector<FT> cands;
  cands.push_back(FT(0));
  for (const Bar& a : A) {
    const FT da = diag_cost(a);
    if (da >= FT(0)) cands.push_back(da);
    for (const Bar& b : B) {
      const FT c = linf_bars(a, b);
      if (c >= FT(0)) cands.push_back(c);
    }
  }
  for (const Bar& b : B) {
    const FT db = diag_cost(b);
    if (db >= FT(0)) cands.push_back(db);
  }
  std::sort(cands.begin(), cands.end());
  cands.erase(std::unique(cands.begin(), cands.end()), cands.end());

  FT ans = FT(0);
  bool found = false;
  int lo = 0, hi = static_cast<int>(cands.size()) - 1;
  while (lo <= hi) {
    const int mid = (lo + hi) / 2;
    if (matching_feasible(A, B, cands[mid])) {
      ans = cands[mid];
      found = true;
      hi = mid - 1;
    } else {
      lo = mid + 1;
    }
  }
  if (!found) {
    if (witness) witness->valid = false;
    return infinite_F();
  }

  if (witness) {
    witness->valid = false;
    auto set_from_grades = [&](const Grade2& p, const Grade2& q, bool same, const FT& s_dummy, const FT& t_dummy) {
      (void)s_dummy;
      (void)t_dummy;
      witness->p = p;
      witness->q = q;
      witness->same_module = same;
      witness->valid = true;
    };
    for (const Bar& a : A) {
      const FT da = diag_cost(a);
      if (da == ans && a.has_death_pt) {
        set_from_grades(a.birth_pt, a.death_pt, true, 0, 0);
        return ans;
      }
    }
    for (const Bar& b : B) {
      const FT db = diag_cost(b);
      if (db == ans && b.has_death_pt) {
        set_from_grades(b.birth_pt, b.death_pt, true, 0, 0);
        return ans;
      }
    }
    for (const Bar& a : A) {
      for (const Bar& b : B) {
        const FT c = linf_bars(a, b);
        if (c != ans) continue;
        if (a.infinite && b.infinite) {
          set_from_grades(a.birth_pt, b.birth_pt, false, 0, 0);
          return ans;
        }
        if (a.infinite || b.infinite) continue;
        if (abs_ft(a.birth - b.birth) == ans) {
          set_from_grades(a.birth_pt, b.birth_pt, false, 0, 0);
          return ans;
        }
        if (a.has_death_pt && b.has_death_pt && abs_ft(a.death - b.death) == ans) {
          set_from_grades(a.death_pt, b.death_pt, false, 0, 0);
          return ans;
        }
      }
    }
  }
  return ans;
}

inline FT evaluate_F(const Presentation& M, const Presentation& N,
                     const FT& s, const FT& t, Witness* witness = nullptr) {
  if (s < FT(0)) {
    if (witness) witness->valid = false;
    return FT(0);
  }
  const std::vector<Bar> dM = slice_diagram(M, s, t);
  const std::vector<Bar> dN = slice_diagram(N, s, t);
  const FT d = bottleneck(dM, dN, witness);
  if (d < FT(0)) return d;
  if (witness && witness->valid) {
    witness->type = delta_type(witness->p, witness->q, s, t);
  }
  return d;
}

inline Line_2 dual_line_of_grade(const Grade2& v) {
  return Line_2(v.x, FT(1), -v.y);
}

inline bool is_degenerate_line(const Line_2& ln) {
  return ln.a() == FT(0) && ln.b() == FT(0);
}

inline void add_line(std::vector<Line_2>& lines, const Line_2& ln) {
  if (is_degenerate_line(ln)) return;
  lines.push_back(ln);
}

inline void add_line(std::vector<X_monotone_curve_2>& curves, const Line_2& ln) {
  if (is_degenerate_line(ln)) return;
  curves.push_back(X_monotone_curve_2(ln));
}

inline bool same_line(const Line_2& a, const Line_2& b) {
  return CGAL::parallel(a, b) && a.has_on(b.point(0));
}

inline void unique_add_line(std::vector<Line_2>& lines, const Line_2& ln) {
  if (is_degenerate_line(ln)) return;
  for (const Line_2& e : lines) {
    if (same_line(e, ln)) return;
  }
  lines.push_back(ln);
}

inline void add_comparison_lines(std::vector<Line_2>& curves,
                                 Grade2 p, Grade2 q, DeltaType tp, const FT& c1,
                                 Grade2 r, Grade2 spt, DeltaType tr, const FT& c2) {
  if (tp == TYPE_IIIB) std::swap(p, q);
  if (tr == TYPE_IIIB) std::swap(r, spt);
  const DeltaType a = (tp == TYPE_IIIB) ? TYPE_IIIA : tp;
  const DeltaType b = (tr == TYPE_IIIB) ? TYPE_IIIA : tr;

  auto add_st = [&](const FT& A, const FT& B, const FT& C) {
    add_line(curves, Line_2(A, B, C));
  };

  if (a == TYPE_I && b == TYPE_I) return;
  if (a == TYPE_II && b == TYPE_II) return;

  if ((a == TYPE_I && b == TYPE_II) || (a == TYPE_II && b == TYPE_I)) {
    const Grade2& Ipq = (a == TYPE_I) ? Grade2{p.y, q.y} : Grade2{r.y, spt.y};
    const Grade2& IIxy = (a == TYPE_II) ? Grade2{p.x, q.x} : Grade2{r.x, spt.x};
    const FT cI = (a == TYPE_I) ? c1 : c2;
    const FT cII = (a == TYPE_II) ? c1 : c2;
    const FT dy = abs_ft(Ipq.x - Ipq.y);
    const FT dx = abs_ft(IIxy.x - IIxy.y);
    if (dx == FT(0)) return;
    const FT s0 = (cI * dy) / (cII * dx);
    if (s0 > FT(0)) add_line(curves, Line_2(Point_2(s0, FT(0)), Point_2(s0, FT(1))));
    return;
  }

  if (a == TYPE_IIIA && b == TYPE_IIIA) {
    add_st(-c1 * q.x + c2 * spt.x, -c1 + c2, c1 * p.y - c2 * r.y);
    add_st(-c1 * q.x - c2 * spt.x, -c1 - c2, c1 * p.y + c2 * r.y);
    return;
  }

  if ((a == TYPE_I && b == TYPE_IIIA) || (a == TYPE_IIIA && b == TYPE_I)) {
    const Grade2 pup = (b == TYPE_IIIA) ? r : p;
    const Grade2 qdn = (b == TYPE_IIIA) ? spt : q;
    const Grade2 gp = (a == TYPE_I) ? p : r;
    const Grade2 gq = (a == TYPE_I) ? q : spt;
    const FT cI = (a == TYPE_I) ? c1 : c2;
    const FT cIII = (b == TYPE_IIIA) ? c2 : c1;
    const FT k = (cI / cIII) * abs_ft(gp.y - gq.y);
    add_st(-qdn.x, FT(-1), pup.y - k);
    add_st(-qdn.x, FT(-1), pup.y + k);
    return;
  }

  if ((a == TYPE_II && b == TYPE_IIIA) || (a == TYPE_IIIA && b == TYPE_II)) {
    const Grade2 pup = (b == TYPE_IIIA) ? r : p;
    const Grade2 qdn = (b == TYPE_IIIA) ? spt : q;
    const Grade2 gp = (a == TYPE_II) ? p : r;
    const Grade2 gq = (a == TYPE_II) ? q : spt;
    const FT cII = (a == TYPE_II) ? c1 : c2;
    const FT cIII = (b == TYPE_IIIA) ? c2 : c1;
    const FT k = (cII / cIII) * abs_ft(gp.x - gq.x);
    add_st(-(qdn.x + k), FT(-1), pup.y);
    add_st(-(qdn.x - k), FT(-1), pup.y);
  }
}

struct LabeledGrade {
  Grade2 g;
  int module;
};

inline FT pair_coeff(const LabeledGrade& a, const LabeledGrade& b) {
  return (a.module == b.module) ? FT(1) / FT(2) : FT(1);
}

inline bool face_interior_point(Arrangement_2::Face_const_handle f, Point_2& out) {
  if (f->is_fictitious()) return false;

  auto try_ccb = [&](Arrangement_2::Ccb_halfedge_const_circulator circ) -> bool {
    Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
    do {
      if (curr->is_fictitious()) {
        ++curr;
        continue;
      }
      const bool src_inf = curr->source()->is_at_open_boundary();
      const bool tgt_inf = curr->target()->is_at_open_boundary();
      Point_2 p, q;
      bool have_two = false;
      if (!src_inf && !tgt_inf) {
        p = curr->source()->point();
        q = curr->target()->point();
        have_two = (p != q);
      } else if (!src_inf || !tgt_inf) {
        const Point_2 base = src_inf ? curr->target()->point() : curr->source()->point();
        const Line_2 ln = curr->curve().supporting_line();
        Point_2 other = ln.point(1);
        if (other == base) other = ln.point(-1);
        p = base;
        q = other;
        have_two = (p != q);
      } else {
        const Line_2 ln = curr->curve().supporting_line();
        p = ln.point(0);
        q = ln.point(1);
        have_two = (p != q);
      }
      if (have_two) {
        const Point_2 mid = CGAL::midpoint(p, q);
        const Vector_2 v = q - p;
        const Vector_2 left(-v.y(), v.x());
        out = mid + left / FT(32);
        if (out.x() > FT(0)) return true;
        out = mid - left / FT(32);
        if (out.x() > FT(0)) return true;
      }
      ++curr;
    } while (curr != circ);
    return false;
  };

  if (f->has_outer_ccb()) {
    for (auto it = f->outer_ccbs_begin(); it != f->outer_ccbs_end(); ++it) {
      if (try_ccb(*it)) return true;
    }
  }
  for (auto it = f->inner_ccbs_begin(); it != f->inner_ccbs_end(); ++it) {
    if (try_ccb(*it)) return true;
  }
  return false;
}

inline std::vector<FT> right_unbounded_slopes(Arrangement_2::Face_const_handle f) {
  std::vector<FT> slopes;
  auto walk = [&](Arrangement_2::Ccb_halfedge_const_circulator circ) {
    Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
    do {
      if (curr->is_fictitious()) {
        ++curr;
        continue;
      }
      const bool src_inf = curr->source()->is_at_open_boundary();
      const bool tgt_inf = curr->target()->is_at_open_boundary();
      if (!src_inf && !tgt_inf) {
        ++curr;
        continue;
      }
      const bool to_plus_inf =
          (tgt_inf && curr->direction() == CGAL::ARR_LEFT_TO_RIGHT) ||
          (src_inf && curr->direction() == CGAL::ARR_RIGHT_TO_LEFT);
      if (!to_plus_inf && !(src_inf && tgt_inf)) {
        ++curr;
        continue;
      }
      const Line_2 ln = curr->curve().supporting_line();
      if (ln.b() == FT(0)) {
        ++curr;
        continue;
      }
      slopes.push_back(-ln.a() / ln.b());
      ++curr;
    } while (curr != circ);
  };
  if (f->has_outer_ccb()) {
    for (auto it = f->outer_ccbs_begin(); it != f->outer_ccbs_end(); ++it) walk(*it);
  }
  std::sort(slopes.begin(), slopes.end());
  slopes.erase(std::unique(slopes.begin(), slopes.end()), slopes.end());
  return slopes;
}

inline void collect_labeled_grades(const Presentation& M, const Presentation& N,
                                   std::vector<LabeledGrade>& grades) {
  grades.clear();
  auto add_unique = [&](const std::vector<Grade2>& pts, int module) {
    std::set<Grade2> seen;
    for (const Grade2& g : pts) {
      if (seen.insert(g).second) {
        LabeledGrade lg;
        lg.g = g;
        lg.module = module;
        grades.push_back(lg);
      }
    }
  };
  add_unique(M.generators, 0);
  add_unique(M.relations, 0);
  add_unique(N.generators, 1);
  add_unique(N.relations, 1);
}

inline void build_arrangement_A(Arrangement_2& arr,
                                const Presentation& M,
                                const Presentation& N,
                                bool /*refine*/) {
  std::vector<LabeledGrade> grades;
  collect_labeled_grades(M, N, grades);

  std::vector<X_monotone_curve_2> curves;
  std::set<Grade2> dual_pts;
  for (const LabeledGrade& lg : grades) {
    if (dual_pts.insert(lg.g).second) {
      add_line(curves, dual_line_of_grade(lg.g));
    }
  }
  curves.push_back(X_monotone_curve_2(Line_2(Point_2(FT(0), FT(0)), Point_2(FT(0), FT(1)))));
  curves.push_back(X_monotone_curve_2(Line_2(Point_2(FT(1), FT(0)), Point_2(FT(1), FT(1)))));
  CGAL::insert(arr, curves.begin(), curves.end());
}

struct PairKey {
  Grade2 p, q;
  bool same_module;
};

inline bool operator<(const PairKey& a, const PairKey& b) {
  if (a.same_module != b.same_module) return a.same_module < b.same_module;
  if (!(a.p == b.p)) return a.p < b.p;
  return a.q < b.q;
}

inline PairKey make_pair_key(Witness w) {
  if (w.q < w.p) std::swap(w.p, w.q);
  PairKey k;
  k.p = w.p;
  k.q = w.q;
  k.same_module = w.same_module;
  return k;
}

inline bool line_intersection_point(const Line_2& a, const Line_2& b, Point_2& out) {
  if (CGAL::parallel(a, b)) return false;
  const auto hit = CGAL::intersection(a, b);
  if (!hit) return false;
  if (const Point_2* p = std::get_if<Point_2>(&*hit)) {
    out = *p;
    return true;
  }
  return false;
}

template<typename PL>
inline bool in_face_closure(const PL& pl, Arrangement_2::Face_const_handle f, const Point_2& pt) {
  const auto obj = pl.locate(pt);
  using FH = Arrangement_2::Face_const_handle;
  using VH = Arrangement_2::Vertex_const_handle;
  using EH = Arrangement_2::Halfedge_const_handle;
  if (const FH* fp = std::get_if<FH>(&obj)) return *fp == f;
  if (const VH* vp = std::get_if<VH>(&obj)) {
    const VH v = *vp;
    if (v->is_isolated()) return v->face() == f;
    auto circ = v->incident_halfedges();
    auto curr = circ;
    do {
      if (curr->face() == f) return true;
      ++curr;
    } while (curr != circ);
    return false;
  }
  if (const EH* ep = std::get_if<EH>(&obj)) {
    return (*ep)->face() == f || (*ep)->twin()->face() == f;
  }
  return false;
}

inline void face_supporting_lines(Arrangement_2::Face_const_handle f, std::vector<Line_2>& lines) {
  auto walk = [&](Arrangement_2::Ccb_halfedge_const_circulator circ) {
    Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
    do {
      if (!curr->is_fictitious()) unique_add_line(lines, curr->curve().supporting_line());
      ++curr;
    } while (curr != circ);
  };
  if (f->has_outer_ccb()) {
    for (auto it = f->outer_ccbs_begin(); it != f->outer_ccbs_end(); ++it) walk(*it);
  }
  for (auto it = f->inner_ccbs_begin(); it != f->inner_ccbs_end(); ++it) walk(*it);
}

inline void comparison_lines_for_witness(std::vector<Line_2>& lines,
                                         const Witness& w,
                                         const std::vector<LabeledGrade>& grades,
                                         const FT& s0, const FT& t0) {
  const FT c1 = w.same_module ? FT(1) / FT(2) : FT(1);
  const int n = static_cast<int>(grades.size());
  for (int i = 0; i < n; ++i) {
    for (int j = i; j < n; ++j) {
      const bool same_pair =
          ((grades[i].g == w.p && grades[j].g == w.q) ||
           (grades[i].g == w.q && grades[j].g == w.p)) &&
          ((grades[i].module == grades[j].module) == w.same_module);
      if (same_pair) continue;
      add_comparison_lines(lines, w.p, w.q, w.type, c1,
                           grades[i].g, grades[j].g,
                           delta_type(grades[i].g, grades[j].g, s0, t0),
                           pair_coeff(grades[i], grades[j]));
    }
  }
  if (w.type == TYPE_IIIA) add_line(lines, Line_2(-w.q.x, FT(-1), w.p.y));
  if (w.type == TYPE_IIIB) add_line(lines, Line_2(-w.p.x, FT(-1), w.q.y));
}

inline void unique_lines(std::vector<Line_2>& lines) {
  std::vector<Line_2> u;
  for (const Line_2& ln : lines) unique_add_line(u, ln);
  lines.swap(u);
}

inline FT maximize_F(const Arrangement_2& arr,
                     const Presentation& M,
                     const Presentation& N,
                     bool refine = true) {
  FT best = FT(0);
  bool any = false;
  bool infinite = false;

  auto consider_finite = [&](const FT& val) {
    if (!any || val > best) {
      best = val;
      any = true;
    }
  };

  // Ω is s > 0. Unmatched infinite bars on s = 0 can appear even when the
  // matching distance is finite, because a horizontal slice never sees
  // grades with larger y. Count infinity only in the open parameter space.
  auto consider_F = [&](const FT& s, const FT& val) {
    if (is_infinite_F(val)) {
      if (s > FT(0)) infinite = true;
      return;
    }
    consider_finite(val);
  };

  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    if (vit->is_at_open_boundary()) continue;
    const Point_2 p = vit->point();
    if (p.x() < FT(0)) continue;
    consider_F(p.x(), evaluate_F(M, N, p.x(), p.y()));
    if (infinite) return infinite_F();
  }

  CGAL::Arr_naive_point_location<Arrangement_2> pl(arr);
  std::vector<LabeledGrade> grades;
  collect_labeled_grades(M, N, grades);

  for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    if (fit->is_fictitious()) continue;
    Point_2 interior;
    if (!face_interior_point(fit, interior)) continue;
    if (interior.x() <= FT(0)) continue;

    Witness w0;
    consider_F(interior.x(), evaluate_F(M, N, interior.x(), interior.y(), &w0));
    if (infinite) return infinite_F();

    if (fit->is_unbounded() && w0.valid &&
        (w0.type == TYPE_IIIA || w0.type == TYPE_IIIB) &&
        interior.x() >= FT(1)) {
      Grade2 p = w0.p, q = w0.q;
      if (w0.type == TYPE_IIIB) std::swap(p, q);
      const FT c = w0.same_module ? FT(1) / FT(2) : FT(1);
      const std::vector<FT> slopes = right_unbounded_slopes(fit);
      for (const FT& r : slopes) consider_finite(c * abs_ft(r + q.x));
    }

    if (!refine || !w0.valid) continue;

    std::set<PairKey> seen;
    std::vector<Witness> queue;
    queue.push_back(w0);
    seen.insert(make_pair_key(w0));

    std::vector<Line_2> boundary;
    face_supporting_lines(fit, boundary);

    for (std::size_t qi = 0; qi < queue.size(); ++qi) {
      const Witness w = queue[qi];
      std::vector<Line_2> lines;
      comparison_lines_for_witness(lines, w, grades, interior.x(), interior.y());
      unique_lines(lines);

      auto consider_pt = [&](const Point_2& pt) {
        if (infinite) return;
        if (pt.x() < FT(0)) return;
        if (!in_face_closure(pl, fit, pt)) return;
        Witness w2;
        consider_F(pt.x(), evaluate_F(M, N, pt.x(), pt.y(), &w2));
        if (infinite || !w2.valid) return;
        const PairKey key = make_pair_key(w2);
        if (seen.insert(key).second) queue.push_back(w2);
      };

      const int nl = static_cast<int>(lines.size());
      for (int i = 0; i < nl; ++i) {
        for (int j = i + 1; j < nl; ++j) {
          Point_2 ip;
          if (line_intersection_point(lines[i], lines[j], ip)) consider_pt(ip);
        }
        for (const Line_2& e : boundary) {
          Point_2 ip;
          if (line_intersection_point(lines[i], e, ip)) consider_pt(ip);
        }
      }
      if (infinite) return infinite_F();
    }
  }

  return any ? best : FT(0);
}

}  // namespace matching
