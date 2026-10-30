#pragma once

// Exact matching distance between two 2-parameter persistence modules,
// following Kerber, Lesnick, Oudot, arXiv:1812.09085, Section 5.
//
// A slice is an affine line ℓ : y = s x + t in the grade plane, with
// parameters (s, t) ∈ Ω = (0, ∞) × ℝ. Restricting a module M to ℓ gives a
// 1-parameter persistence module M^ℓ and a barcode. The matching distance is
//
//     d_match(M, N) = sup_ℓ  F(ℓ),    F(ℓ) = w(ℓ) · d_B(M^ℓ, N^ℓ).
//
// Weighted grades (below) fold the slice weight w(ℓ) into the 1-parameter
// filtration, so the ordinary bottleneck of those diagrams is already F.
//
// The parameter plane is subdivided by the coarse arrangement A0: the dual
// line of every input grade, plus the verticals s = 0 and s = 1. On each
// cell of A0 the combinatorial type of F is constant except across
// comparison lines of the realizing pair.
//
// maximize_F samples F at A0 vertices, a true interior point of each face,
// s → +∞ limits, and (if refine) finite-edge midpoints plus a few extra
// interior points, giving a lower bound δ. Each bounded A0 face then gets
// the Kerber–Nigmetov variation bound
//     U(f) = F(center) + v(M, f) + v(N, f),
// with v the max change of a weighted push over the vertices of f.
// Extrema of wpush sit at vertices: ∂/∂t of (y−t)/s is never zero, and the
// restriction to a straight edge is monotonic. Unbounded faces have
// infinite diameter and stay red. Green faces (U(f) ≤ δ) are skipped.
//
// On each remaining red face, comparison lines of the current realizing
// pair (O(n²), never the global O(n⁴) overlay) are shuffled and inserted
// incrementally into a local arrangement of that face. New vertices update
// δ; split subfaces with U ≤ δ are marked green; green–green edges are
// deleted. `--no-refine` stops after the A0 sample.

#include <algorithm>
#include <functional>
#include <iostream>
#include <random>
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
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/intersections.h>
#include <CGAL/number_utils.h>

namespace matching {

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef CGAL::Arr_linear_traits_2<Kernel> Traits_2;
typedef Traits_2::Point_2 Point_2;   // (s, t) in slice-parameter space
typedef Traits_2::Line_2 Line_2;     // a s + b t + c = 0
typedef Traits_2::Ray_2 Ray_2;
typedef Traits_2::Segment_2 Segment_2;
typedef Traits_2::Vector_2 Vector_2;
typedef Traits_2::X_monotone_curve_2 X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
// Local overlay of one A0 face: face data is FACE_RED or FACE_GREEN.
typedef CGAL::Arr_face_extended_dcel<Traits_2, int> Face_colored_dcel;
typedef CGAL::Arrangement_2<Traits_2, Face_colored_dcel> Colored_arrangement;
const int FACE_RED = 0;
const int FACE_GREEN = 1;

// a0_lines: duals plus s=0 and s=1 inserted into the global A0.
// lines_introduced: a0_lines plus comparison lines that actually entered
// a local overlay (one count per geometric line, not per clipped piece).
// Green-green deletions are remove_edge calls that merge two green faces.
struct ArrangementStats {
  std::size_t a0_lines = 0;
  std::size_t lines_introduced = 0;
  std::size_t green_green_deletions = 0;
};

// Bigrade in the module's (x, y) plane.
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

// Finite presentation: generators, relations, and one column per relation
// listing the generators it hits (GF(2), 0-based indices).
struct Presentation {
  std::vector<Grade2> generators;
  std::vector<Grade2> relations;
  std::vector<std::vector<int> > columns;
};

// Direct sum of presentations (AIDA indecomposables of one module).
// Matching distance does not factor over summands; F uses the union of
// their slice barcodes. A0 uses the union of their distinct grades.
struct Module {
  std::vector<Presentation> summands;
};

inline Module as_module(Presentation P) {
  Module M;
  M.summands.push_back(std::move(P));
  return M;
}

inline FT abs_ft(const FT& v) { return v < 0 ? -v : v; }

inline FT max_ft(const FT& a, const FT& b) { return a < b ? b : a; }

// Exact FT has no +∞. Negative bottleneck / F means an infinite value
// (slice diagrams with different numbers of infinite bars).
inline FT infinite_F() { return FT(-1); }
inline bool is_infinite_F(const FT& d) { return d < FT(0); }

// How a realizing pair of grades (p, q) sits relative to the slice y = s x + t.
// The formula for their contribution δ_{p,q} to F depends on this type
// (paper, §3–4): both above (I), both below (II), or one on each side (III).
enum DeltaType { TYPE_I, TYPE_II, TYPE_IIIA, TYPE_IIIB };

inline bool point_on_or_above_slice(const Grade2& p, const FT& s, const FT& t) {
  return p.y >= s * p.x + t;
}

inline DeltaType delta_type(const Grade2& p, const Grade2& q, const FT& s, const FT& t) {
  const bool p_above = point_on_or_above_slice(p, s, t);
  const bool q_above = point_on_or_above_slice(q, s, t);
  if (p_above && q_above) return TYPE_I;
  if (!p_above && !q_above) return TYPE_II;
  if (p_above) return TYPE_IIIA;  // p above, q below
  return TYPE_IIIB;               // p below, q above
}

// Push of a 2-parameter grade onto the slice, times the slice weight w(ℓ).
// Bottleneck of diagrams built from these values is F itself, so the
// arrangement never has to carry w(ℓ) separately. The two formulas switch
// at s = 1, where the weight expression in the paper changes.
inline FT weighted_grade(const Grade2& p, const FT& s, const FT& t) {
  const bool above = point_on_or_above_slice(p, s, t);
  if (s >= FT(1)) {
    return above ? (p.y - t) / s : p.x;
  }
  return above ? (p.y - t) : s * p.x;
}

// One bar of a slice barcode, with the original 2-parameter grades that
// produced the birth (and death, if finite). Those grades are the witness
// pair used later to build comparison lines.
struct Bar {
  FT birth;
  FT death;
  bool infinite;
  Grade2 birth_pt;
  Grade2 death_pt;
  bool has_death_pt;
};

// Persistence of a presentation whose generators/relations have already been
// assigned 1-parameter grades (weighted slice grades, or their s → +∞ limits).
inline std::vector<Bar> slice_diagram_rank1(const Presentation& P,
                                            const std::vector<FT>& gen_val,
                                            const std::vector<FT>& rel_val) {
  std::vector<Bar> bars;
  Bar b;
  b.birth = gen_val[0];
  b.birth_pt = P.generators[0];
  bool kills = false;
  if (!P.relations.empty() && !P.columns.empty()) {
    for (int g : P.columns[0]) {
      if (g == 0) {
        kills = true;
        break;
      }
    }
  }
  if (!kills) {
    b.infinite = true;
    b.death = FT(0);
    b.has_death_pt = false;
    bars.push_back(b);
    return bars;
  }
  b.infinite = false;
  b.death = rel_val[0];
  b.death_pt = P.relations[0];
  b.has_death_pt = true;
  if (b.birth < b.death) bars.push_back(b);
  return bars;
}

inline std::vector<Bar> slice_diagram_from_values(const Presentation& P,
                                                  const std::vector<FT>& gen_val,
                                                  const std::vector<FT>& rel_val) {
  const int n_gen = static_cast<int>(P.generators.size());
  const int n_rel = static_cast<int>(P.relations.size());
  std::vector<Bar> bars;
  if (n_gen == 0) return bars;
  if (n_gen == 1 && n_rel <= 1) return slice_diagram_rank1(P, gen_val, rel_val);

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

  // PHAT matrix in filtration order: columns 0 .. n_gen-1 are generators
  // (empty boundary, dim 0); columns n_gen .. are relations (dim 1), with
  // faces rewritten to the sorted generator indices.
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
  // Finite bars: a generator paired with a later relation.
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
  // Unpaired generators are infinite bars (classes that never die on ℓ).
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

// Persistence of the presentation along the slice (s, t): sort generators
// and relations by weighted grade, reduce the resulting 1-parameter complex
// with PHAT, and read off finite and infinite bars.
inline std::vector<Bar> slice_diagram(const Presentation& P, const FT& s, const FT& t) {
  const int n_gen = static_cast<int>(P.generators.size());
  const int n_rel = static_cast<int>(P.relations.size());
  std::vector<FT> gen_val(n_gen), rel_val(n_rel);
  for (int i = 0; i < n_gen; ++i) gen_val[i] = weighted_grade(P.generators[i], s, t);
  for (int i = 0; i < n_rel; ++i) rel_val[i] = weighted_grade(P.relations[i], s, t);
  return slice_diagram_from_values(P, gen_val, rel_val);
}

inline std::vector<Bar> slice_diagram(const Module& M, const FT& s, const FT& t) {
  std::vector<Bar> bars;
  for (const Presentation& P : M.summands) {
    std::vector<Bar> part = slice_diagram(P, s, t);
    bars.insert(bars.end(), part.begin(), part.end());
  }
  return bars;
}

// Realizing pair of F: the two grades whose weighted distance equals the
// bottleneck. same_module is true when the pair is a bar matched to the
// diagonal (birth and death of one module) rather than a match across M and N.
struct Witness {
  Grade2 p, q;
  bool same_module;
  DeltaType type;
  bool valid;
};

// L∞ cost of matching two bars. Negative means the pair is forbidden
// (an infinite bar cannot match a finite one).
inline FT linf_bars(const Bar& a, const Bar& b) {
  if (a.infinite && b.infinite) return abs_ft(a.birth - b.birth);
  if (a.infinite || b.infinite) return FT(-1); // infeasible
  return max_ft(abs_ft(a.birth - b.birth), abs_ft(a.death - b.death));
}

inline FT diag_cost(const Bar& a) {
  if (a.infinite) return FT(-1);  // infinite bars cannot go to the diagonal
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

// True iff there is a partial matching of cost ≤ eps: every bar that is
// longer than 2 eps (so cannot go to the diagonal) is matched to the other
// diagram. DFS augmenting paths on the bipartite graph of allowed pairs.
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

// Bottleneck distance of two slice diagrams. Candidate costs are the
// finitely many pairwise L∞ values and diagonal costs; binary search the
// smallest feasible one. If witness is non-null, fill it with a pair of
// grades that attains that cost (used to spawn comparison lines).
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

  // Recover a realizing pair: prefer a diagonal match (birth/death of one
  // bar), else a cross-match of two bars whose L∞ cost equals the bottleneck.
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

// F(s, t) = bottleneck of the two weighted slice diagrams. Negative slope
// is outside Ω and contributes 0. A negative return is the infinity sentinel.
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

inline FT evaluate_F(const Module& M, const Module& N,
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

// Weighted grade as s → +∞ along a ray with t/s → r, on the A0 face whose
// above/below partition is that of (s_ref, t_ref). Points below the slice
// keep their x-coordinate; points above all collapse to −r.
inline FT grade_at_s_infinity(const Grade2& p, const FT& r,
                              const FT& s_ref, const FT& t_ref) {
  return point_on_or_above_slice(p, s_ref, t_ref) ? -r : p.x;
}

// Bottleneck of the two slice diagrams in that s → +∞ limit. This is the
// actual lim F along the ray, not the closed-form c|r + q_x| of an interior
// witness pair — that pair need not realize F out to infinity on an A0 face.
inline FT evaluate_F_s_infinity(const Presentation& M, const Presentation& N,
                                const FT& r, const FT& s_ref, const FT& t_ref,
                                Witness* witness = nullptr) {
  auto diagram = [&](const Presentation& P) {
    const int n_gen = static_cast<int>(P.generators.size());
    const int n_rel = static_cast<int>(P.relations.size());
    std::vector<FT> gen_val(n_gen), rel_val(n_rel);
    for (int i = 0; i < n_gen; ++i)
      gen_val[i] = grade_at_s_infinity(P.generators[i], r, s_ref, t_ref);
    for (int i = 0; i < n_rel; ++i)
      rel_val[i] = grade_at_s_infinity(P.relations[i], r, s_ref, t_ref);
    return slice_diagram_from_values(P, gen_val, rel_val);
  };
  const FT d = bottleneck(diagram(M), diagram(N), witness);
  if (d < FT(0)) return d;
  if (witness && witness->valid) {
    // The limiting slice is vertical of slope s = ∞; treat as s ≥ 1.
    witness->type = delta_type(witness->p, witness->q, s_ref, t_ref);
  }
  return d;
}

inline std::vector<Bar> slice_diagram_s_infinity(const Presentation& P, const FT& r,
                                                 const FT& s_ref, const FT& t_ref) {
  const int n_gen = static_cast<int>(P.generators.size());
  const int n_rel = static_cast<int>(P.relations.size());
  std::vector<FT> gen_val(n_gen), rel_val(n_rel);
  for (int i = 0; i < n_gen; ++i)
    gen_val[i] = grade_at_s_infinity(P.generators[i], r, s_ref, t_ref);
  for (int i = 0; i < n_rel; ++i)
    rel_val[i] = grade_at_s_infinity(P.relations[i], r, s_ref, t_ref);
  return slice_diagram_from_values(P, gen_val, rel_val);
}

inline FT evaluate_F_s_infinity(const Module& M, const Module& N,
                                const FT& r, const FT& s_ref, const FT& t_ref,
                                Witness* witness = nullptr) {
  auto diagram = [&](const Module& Mod) {
    std::vector<Bar> bars;
    for (const Presentation& P : Mod.summands) {
      std::vector<Bar> part = slice_diagram_s_infinity(P, r, s_ref, t_ref);
      bars.insert(bars.end(), part.begin(), part.end());
    }
    return bars;
  };
  const FT d = bottleneck(diagram(M), diagram(N), witness);
  if (d < FT(0)) return d;
  if (witness && witness->valid) {
    witness->type = delta_type(witness->p, witness->q, s_ref, t_ref);
  }
  return d;
}

// Dual of a grade v = (v_x, v_y): the line t = v_y − s v_x in parameter
// space. Crossing it changes whether v lies above or below the slice.
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

// Insert ln unless a coincident line is already present.
inline void unique_add_line(std::vector<Line_2>& lines, const Line_2& ln) {
  if (is_degenerate_line(ln)) return;
  for (const Line_2& e : lines) {
    if (same_line(e, ln)) return;
  }
  lines.push_back(ln);
}

// Locus in (s, t) where two candidate pairs have equal δ, i.e. where the
// realizing pair of F can switch. c1, c2 are 1 (match across modules) or
// 1/2 (match to the diagonal). Type IIIB is rewritten as IIIA by swapping
// the two grades. Same-type I/I and II/II never cross: those δ are constant.
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

  // Type I vs II: δ_I is a vertical gap, δ_II a horizontal gap times s.
  // They meet on a vertical line s = s0.
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
    // Type I vs III: δ_I is constant, δ_III is linear; they meet on two
    // translates of the type-III dual line.
    const FT k = (cI / cIII) * abs_ft(gp.y - gq.y);
    add_st(-qdn.x, FT(-1), pup.y - k);
    add_st(-qdn.x, FT(-1), pup.y + k);
    return;
  }

  // Type II vs III: same idea, k is a horizontal gap times s.
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
  int module;  // 0 = M, 1 = N
};

// 1/2 when both grades belong to the same module (diagonal match).
inline FT pair_coeff(const LabeledGrade& a, const LabeledGrade& b) {
  return (a.module == b.module) ? FT(1) / FT(2) : FT(1);
}

template<typename PL, typename FaceHandle>
inline bool in_face_interior(const PL& pl, FaceHandle f, const Point_2& pt) {
  const auto obj = pl.locate(pt);
  // locate() stores const handles; a mutable Face_handle is a different type.
  using FH = typename PL::Face_const_handle;
  if (const FH* fp = std::get_if<FH>(&obj)) return *fp == f;
  return false;
}

inline Vector_2 l1_unit(const Vector_2& v) {
  const FT n = FT(1) + abs_ft(v.x()) + abs_ft(v.y());
  return Vector_2(v.x() / n, v.y() / n);
}

// Direction of h from source to target, using the supporting line (works
// when one or both ends sit at infinity).
template<typename HalfedgeHandle>
inline Vector_2 oriented_halfedge_dir(HalfedgeHandle h) {
  Vector_2 dir = h->curve().supporting_line().to_vector();
  if (dir.x() < FT(0) || (dir.x() == FT(0) && dir.y() < FT(0)))
    dir = -dir;
  if (h->direction() == CGAL::ARR_RIGHT_TO_LEFT)
    dir = -dir;
  return dir;
}

// A point strictly inside the face with s > 0. Step a short L1-normalized
// offset to the left of a supporting edge (the face is on the left) and
// keep the sample only if point location returns this open face. Offsets
// proportional to ln.point(1) on unbounded edges were large enough to
// land in a neighbouring cell; evaluate_F_s_infinity then mixed that
// cell's above/below split with this face's outer slopes and overestimated.
template<typename PL, typename FaceHandle>
inline bool face_interior_point(const PL& pl, FaceHandle f, Point_2& out) {
  if (f->is_fictitious()) return false;

  static const int k_off[] = {8, 32, 128, 512, 2048, 8192, 32768, 131072, 524288};

  auto try_from = [&](const Point_2& base, const Vector_2& left) -> bool {
    const Vector_2 u = l1_unit(left);
    if (u.x() == FT(0) && u.y() == FT(0)) return false;
    for (int k : k_off) {
      const Point_2 p = base + u / FT(k);
      if (p.x() <= FT(0)) continue;
      if (in_face_interior(pl, f, p)) {
        out = p;
        return true;
      }
    }
    return false;
  };

  auto try_ccb = [&](auto circ) -> bool {
    auto curr = circ;
    do {
      if (curr->is_fictitious()) {
        ++curr;
        continue;
      }
      const Vector_2 dir = oriented_halfedge_dir(curr);
      const Vector_2 left(-dir.y(), dir.x());
      const bool src_inf = curr->source()->is_at_open_boundary();
      const bool tgt_inf = curr->target()->is_at_open_boundary();

      if (!src_inf && !tgt_inf) {
        const Point_2 p = curr->source()->point();
        const Point_2 q = curr->target()->point();
        if (p != q && try_from(CGAL::midpoint(p, q), left)) return true;
      } else if (!src_inf || !tgt_inf) {
        const Point_2 base = src_inf ? curr->target()->point() : curr->source()->point();
        const Vector_2 along = l1_unit(src_inf ? -dir : dir);
        for (int k : k_off) {
          if (try_from(base + along / FT(k), left)) return true;
        }
      } else {
        const Line_2 ln = curr->curve().supporting_line();
        Point_2 on_line;
        bool have = false;
        if (ln.b() != FT(0)) {
          const FT s = FT(1);
          on_line = Point_2(s, -(ln.a() * s + ln.c()) / ln.b());
          have = true;
        } else if (ln.a() != FT(0)) {
          const FT s = -ln.c() / ln.a();
          if (s > FT(0)) {
            on_line = Point_2(s, FT(0));
            have = true;
          }
        }
        if (have && try_from(on_line, left)) return true;
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

// Slopes s of edges of an unbounded face that run to s → +∞. Type-III
// δ grows linearly in s on those rays, and the paper evaluates the limit
// along each such outer direction (Lemma 10).
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

// Distinct generator and relation grades of both modules, tagged with
// which module they came from (needed for the 1 vs 1/2 coefficient).
inline void collect_labeled_grades(const Module& M, const Module& N,
                                   std::vector<LabeledGrade>& grades) {
  grades.clear();
  auto add_module = [&](const Module& Mod, int module) {
    std::set<Grade2> seen;
    auto add = [&](const std::vector<Grade2>& pts) {
      for (const Grade2& g : pts) {
        if (seen.insert(g).second) {
          LabeledGrade lg;
          lg.g = g;
          lg.module = module;
          grades.push_back(lg);
        }
      }
    };
    for (const Presentation& P : Mod.summands) {
      add(P.generators);
      add(P.relations);
    }
  };
  add_module(M, 0);
  add_module(N, 1);
}

inline void collect_labeled_grades(const Presentation& M, const Presentation& N,
                                   std::vector<LabeledGrade>& grades) {
  collect_labeled_grades(as_module(M), as_module(N), grades);
}

// Coarse arrangement A0: one dual line per distinct input grade, plus
// s = 0 and s = 1. The unused refine flag is kept so the call site can
// pass the same boolean it later gives to maximize_F; comparison lines
// are not inserted here.
inline void build_arrangement_A(Arrangement_2& arr,
                                const Module& M,
                                const Module& N,
                                bool /*refine*/,
                                ArrangementStats* stats = nullptr) {
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
  if (stats) {
    stats->a0_lines += curves.size();
    stats->lines_introduced += curves.size();
  }
  CGAL::insert(arr, curves.begin(), curves.end());
}

inline void build_arrangement_A(Arrangement_2& arr,
                                const Presentation& M,
                                const Presentation& N,
                                bool refine,
                                ArrangementStats* stats = nullptr) {
  build_arrangement_A(arr, as_module(M), as_module(N), refine, stats);
}

// Dedup key for realizing pairs while walking comparison vertices of a face.
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

// True if pt lies in the closed face (interior, a bounding vertex, or a
// bounding edge). Comparison-line intersections on the boundary still
// belong to the cell for maximization.
template<typename PL, typename FaceHandle>
inline bool in_face_closure(const PL& pl, FaceHandle f, const Point_2& pt) {
  const auto obj = pl.locate(pt);
  using FH = typename PL::Face_const_handle;
  using VH = typename PL::Vertex_const_handle;
  using EH = typename PL::Halfedge_const_handle;
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

// Supporting lines of every non-fictitious edge of the face (outer and holes).
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

// Comparison lines of the current realizing pair (w.p, w.q) against every
// other grade pair, plus the type-III supporting line of w itself. O(n²)
// lines, local to this witness and this face.
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
  // Dual of the type-III pair: the line through which p and q lie on
  // opposite sides of the slice. Needed so the max on the cell can sit
  // where that pair's own δ has a kink.
  if (w.type == TYPE_IIIA) add_line(lines, Line_2(-w.q.x, FT(-1), w.p.y));
  if (w.type == TYPE_IIIB) add_line(lines, Line_2(-w.p.x, FT(-1), w.q.y));
}

inline void unique_lines(std::vector<Line_2>& lines) {
  std::vector<Line_2> u;
  for (const Line_2& ln : lines) unique_add_line(u, ln);
  lines.swap(u);
}

template<typename FaceHandle, typename Fn>
inline void for_each_face_halfedge(FaceHandle f, Fn&& fn) {
  auto walk = [&](auto circ) {
    auto curr = circ;
    do {
      fn(curr);
      ++curr;
    } while (curr != circ);
  };
  if (f->has_outer_ccb()) {
    for (auto it = f->outer_ccbs_begin(); it != f->outer_ccbs_end(); ++it) walk(*it);
  }
  for (auto it = f->inner_ccbs_begin(); it != f->inner_ccbs_end(); ++it) walk(*it);
}

template<typename FaceHandle>
inline void face_finite_vertices(FaceHandle f, std::vector<Point_2>& verts) {
  verts.clear();
  for_each_face_halfedge(f, [&](auto h) {
    if (h->source()->is_at_open_boundary()) return;
    const Point_2& p = h->source()->point();
    if (p.x() < FT(0)) return;
    verts.push_back(p);
  });
}

template<typename FaceHandle>
inline void face_boundary_curves(FaceHandle f, std::vector<X_monotone_curve_2>& curves) {
  curves.clear();
  for_each_face_halfedge(f, [&](auto h) {
    if (!h->is_fictitious()) curves.push_back(h->curve());
  });
}

// Max change of a weighted push of any grade of Mod over the given vertices.
// On an A0 cell the above/below partition is constant and s = 1 is a wall, so
// each wpush is (y−t)/s, y−t, s x, or constant. None of those has an interior
// critical point on a convex polygon, and each is monotonic on a straight
// edge, so vertex values give the exact variation (KN L-bound, Lemma 2).
inline FT module_variation(const Module& Mod, const std::vector<Point_2>& verts) {
  if (verts.empty()) return FT(0);
  FT vmax = FT(0);
  auto acc = [&](const Grade2& g) {
    FT lo = weighted_grade(g, verts[0].x(), verts[0].y());
    FT hi = lo;
    for (std::size_t i = 1; i < verts.size(); ++i) {
      const FT w = weighted_grade(g, verts[i].x(), verts[i].y());
      if (w < lo) lo = w;
      if (w > hi) hi = w;
    }
    const FT v = hi - lo;
    if (v > vmax) vmax = v;
  };
  for (const Presentation& P : Mod.summands) {
    for (const Grade2& g : P.generators) acc(g);
    for (const Grade2& g : P.relations) acc(g);
  }
  return vmax;
}

inline bool kn_face_upper_bound(const Module& M, const Module& N,
                                const FT& F_center,
                                const std::vector<Point_2>& verts,
                                FT& ub) {
  if (is_infinite_F(F_center) || verts.size() < 2) return false;
  ub = F_center + module_variation(M, verts) + module_variation(N, verts);
  return true;
}

inline bool face_centroid(const std::vector<Point_2>& verts, Point_2& out) {
  if (verts.size() < 3) return false;
  FT sx = FT(0), sy = FT(0);
  for (const Point_2& p : verts) {
    sx += p.x();
    sy += p.y();
  }
  const FT n = FT(static_cast<int>(verts.size()));
  out = Point_2(sx / n, sy / n);
  return out.x() > FT(0);
}

template<typename HalfedgeHandle>
inline bool point_on_halfedge(HalfedgeHandle h, const Point_2& ip) {
  const X_monotone_curve_2& cv = h->curve();
  if (cv.is_segment()) return cv.segment().has_on(ip);
  if (cv.is_ray()) return cv.ray().has_on(ip);
  if (cv.is_line()) return cv.line().has_on(ip);
  return false;
}

template<typename HalfedgeHandle>
inline bool point_on_or_left(HalfedgeHandle h, const Point_2& pt) {
  const Vector_2 dir = oriented_halfedge_dir(h);
  Point_2 base;
  if (!h->source()->is_at_open_boundary()) base = h->source()->point();
  else if (!h->target()->is_at_open_boundary()) base = h->target()->point();
  else base = h->curve().supporting_line().point(0);
  const Vector_2 w = pt - base;
  return dir.x() * w.y() - dir.y() * w.x() >= FT(0);
}

// A0 cells and red subcells of a line overlay are convex, with the face on
// the left of every bounding halfedge.
template<typename FaceHandle>
inline bool point_in_convex_face(FaceHandle f, const Point_2& pt) {
  bool any = false;
  bool ok = true;
  for_each_face_halfedge(f, [&](auto h) {
    if (h->is_fictitious()) return;
    any = true;
    if (!point_on_or_left(h, pt)) ok = false;
  });
  return any && ok;
}

inline FT line_param(const Line_2& L, const Point_2& p) {
  const Vector_2 d = L.to_vector();
  const Vector_2 w = p - L.point(0);
  const FT den = d.x() * d.x() + d.y() * d.y();
  if (den == FT(0)) return FT(0);
  return (d.x() * w.x() + d.y() * w.y()) / den;
}

// L ∩ f as one segment, ray, or line. Convex faces give at most one piece.
template<typename FaceHandle>
inline void clip_line_to_face(FaceHandle f, const Line_2& L,
                              std::vector<X_monotone_curve_2>& pieces) {
  if (is_degenerate_line(L)) return;
  std::vector<std::pair<FT, Point_2> > hits;
  bool coincident = false;
  for_each_face_halfedge(f, [&](auto h) {
    if (h->is_fictitious() || coincident) return;
    const Line_2 e = h->curve().supporting_line();
    if (CGAL::parallel(L, e) && L.has_on(e.point(0))) {
      coincident = true;
      return;
    }
    Point_2 ip;
    if (!line_intersection_point(L, e, ip)) return;
    if (!point_on_halfedge(h, ip)) return;
    hits.push_back(std::make_pair(line_param(L, ip), ip));
  });
  if (coincident) return;

  std::sort(hits.begin(), hits.end(),
            [](const std::pair<FT, Point_2>& a, const std::pair<FT, Point_2>& b) {
              return a.first < b.first;
            });
  std::vector<std::pair<FT, Point_2> > uniq;
  for (const auto& h : hits) {
    if (!uniq.empty() && uniq.back().second == h.second) continue;
    uniq.push_back(h);
  }
  hits.swap(uniq);

  if (hits.size() >= 2) {
    const Point_2& a = hits.front().second;
    const Point_2& b = hits.back().second;
    if (a != b) pieces.push_back(X_monotone_curve_2(Segment_2(a, b)));
    return;
  }

  const Vector_2 dir = L.to_vector();
  const Vector_2 u = l1_unit(dir);
  if (u.x() == FT(0) && u.y() == FT(0)) return;

  auto step_inside = [&](const Point_2& base, const Vector_2& step) -> bool {
    static const int k_off[] = {8, 32, 128, 512, 2048};
    for (int k : k_off) {
      if (point_in_convex_face(f, base + step / FT(k))) return true;
    }
    return false;
  };

  if (hits.size() == 1) {
    const Point_2& p = hits[0].second;
    const bool plus = step_inside(p, u);
    const bool minus = step_inside(p, -u);
    if (plus && minus) pieces.push_back(X_monotone_curve_2(L));
    else if (plus) pieces.push_back(X_monotone_curve_2(Ray_2(p, dir)));
    else if (minus) pieces.push_back(X_monotone_curve_2(Ray_2(p, -dir)));
    return;
  }

  Point_2 test;
  bool have = false;
  if (L.b() != FT(0)) {
    test = Point_2(FT(1), -(L.a() + L.c()) / L.b());
    have = true;
  } else if (L.a() != FT(0)) {
    const FT s = -L.c() / L.a();
    if (s > FT(0)) {
      test = Point_2(s, FT(0));
      have = true;
    }
  }
  if (have && point_in_convex_face(f, test))
    pieces.push_back(X_monotone_curve_2(L));
}

inline bool local_has_red(const Colored_arrangement& local) {
  for (auto fit = local.faces_begin(); fit != local.faces_end(); ++fit) {
    if (fit->is_fictitious()) continue;
    if (fit->data() == FACE_RED) return true;
  }
  return false;
}

inline void merge_green_green(Colored_arrangement& local, ArrangementStats& stats) {
  for (;;) {
    bool killed = false;
    for (auto eit = local.edges_begin(); eit != local.edges_end(); ++eit) {
      if (eit->face()->is_fictitious() || eit->twin()->face()->is_fictitious())
        continue;
      if (eit->face()->data() == FACE_GREEN &&
          eit->twin()->face()->data() == FACE_GREEN) {
        local.remove_edge(eit);
        ++stats.green_green_deletions;
        killed = true;
        break;
      }
    }
    if (!killed) return;
  }
}

inline bool init_local_from_a0_face(Colored_arrangement& local,
                                    Arrangement_2::Face_const_handle a0,
                                    const Point_2& interior) {
  std::vector<X_monotone_curve_2> curves;
  face_boundary_curves(a0, curves);
  if (curves.empty()) return false;
  CGAL::insert(local, curves.begin(), curves.end());
  CGAL::Arr_naive_point_location<Colored_arrangement> lpl(local);
  Colored_arrangement::Face_handle red = local.faces_begin();
  bool found = false;
  for (auto fit = local.faces_begin(); fit != local.faces_end(); ++fit) {
    if (fit->is_fictitious()) continue;
    fit->set_data(FACE_GREEN);
    if (in_face_interior(lpl, fit, interior)) {
      red = fit;
      found = true;
    }
  }
  if (!found) return false;
  red->set_data(FACE_RED);
  return true;
}

template<typename Consider>
inline void recolor_local_faces(Colored_arrangement& local,
                                const Module& M, const Module& N,
                                const FT& best,
                                Consider& consider_F,
                                bool& infinite) {
  CGAL::Arr_naive_point_location<Colored_arrangement> lpl(local);
  for (auto fit = local.faces_begin(); fit != local.faces_end(); ++fit) {
    if (fit->is_fictitious() || fit->data() == FACE_GREEN) continue;
    // Unbounded subfaces can have infinite wpush variation (vertical
    // t → ±∞ rays); they stay red until a later split bounds them.
    if (fit->is_unbounded()) continue;
    Point_2 ip;
    if (!face_interior_point(lpl, fit, ip)) continue;
    if (ip.x() <= FT(0)) {
      fit->set_data(FACE_GREEN);
      continue;
    }
    Witness w;
    const FT fc = evaluate_F(M, N, ip.x(), ip.y(), &w);
    consider_F(ip.x(), fc);
    if (infinite) return;
    if (is_infinite_F(fc)) continue;
    std::vector<Point_2> verts;
    face_finite_vertices(fit, verts);
    FT ub;
    if (kn_face_upper_bound(M, N, fc, verts, ub) && !(best < ub))
      fit->set_data(FACE_GREEN);
  }
}

// Previous all-pairs intersection of a witness's comparison lines in one
// A0 face. Used if the local overlay of that face cannot be built.
template<typename Consider>
inline void refine_face_pairwise(
    Arrangement_2::Face_const_handle f,
    const Point_2& interior,
    const std::vector<Witness>& seeds,
    const Module& M, const Module& N,
    const std::vector<LabeledGrade>& grades,
    const CGAL::Arr_naive_point_location<Arrangement_2>& pl,
    Consider& consider_F,
    bool& infinite) {
  std::set<PairKey> seen;
  std::vector<Witness> queue;
  for (const Witness& w : seeds) {
    if (!w.valid) continue;
    if (seen.insert(make_pair_key(w)).second) queue.push_back(w);
  }
  if (queue.empty()) return;

  std::vector<Line_2> boundary;
  face_supporting_lines(f, boundary);

  for (std::size_t qi = 0; qi < queue.size(); ++qi) {
    const Witness w = queue[qi];
    std::vector<Line_2> lines;
    comparison_lines_for_witness(lines, w, grades, interior.x(), interior.y());
    unique_lines(lines);

    auto consider_pt = [&](const Point_2& pt) {
      if (infinite) return;
      if (pt.x() < FT(0)) return;
      if (!in_face_closure(pl, f, pt)) return;
      Witness w2;
      consider_F(pt.x(), evaluate_F(M, N, pt.x(), pt.y(), &w2));
      if (infinite || !w2.valid) return;
      if (seen.insert(make_pair_key(w2)).second) queue.push_back(w2);
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
    if (infinite) return;
  }
}

// Shuffle the current witness's comparison lines and insert them one at a
// time into a local overlay of this A0 face. Skip a line whose zone is
// entirely green. Boundary hits of a red split raise δ; new subfaces with
// U ≤ δ turn green; green–green edges are deleted. New realizing pairs
// enqueue their own O(n²) lines. Never materializes the global O(n⁴) list.
template<typename Consider>
inline void refine_red_a0_face(
    Arrangement_2::Face_const_handle a0_face,
    const Point_2& interior,
    const std::vector<Witness>& seeds,
    const Module& M, const Module& N,
    const std::vector<LabeledGrade>& grades,
    const CGAL::Arr_naive_point_location<Arrangement_2>& pl,
    bool has_a0_ub, const FT& a0_ub,
    Consider& consider_F,
    FT& best,
    bool& infinite,
    ArrangementStats& stats) {
  Colored_arrangement local;
  if (!init_local_from_a0_face(local, a0_face, interior)) {
    refine_face_pairwise(a0_face, interior, seeds, M, N, grades, pl,
                         consider_F, infinite);
    return;
  }

  std::set<PairKey> seen;
  std::vector<Line_2> pending;
  std::vector<Line_2> used;
  std::mt19937 rng(1);

  auto enqueue_witness_lines = [&](const Witness& w) {
    if (!w.valid) return;
    std::vector<Line_2> lines;
    comparison_lines_for_witness(lines, w, grades, interior.x(), interior.y());
    unique_lines(lines);
    std::shuffle(lines.begin(), lines.end(), rng);
    for (const Line_2& L : lines) unique_add_line(pending, L);
  };

  for (const Witness& w : seeds) {
    if (!w.valid) continue;
    if (seen.insert(make_pair_key(w)).second) enqueue_witness_lines(w);
  }

  auto a0_already_green = [&]() {
    return has_a0_ub && !(best < a0_ub);
  };

  auto eval_hit = [&](const Point_2& pt) {
    if (infinite || pt.x() < FT(0)) return;
    Witness w2;
    consider_F(pt.x(), evaluate_F(M, N, pt.x(), pt.y(), &w2));
    if (infinite || !w2.valid) return;
    if (seen.insert(make_pair_key(w2)).second) enqueue_witness_lines(w2);
  };

  while (!pending.empty() && !infinite) {
    if (a0_already_green()) break;
    if (!local_has_red(local)) break;

    const Line_2 L = pending.back();
    pending.pop_back();
    bool already = false;
    for (const Line_2& e : used) {
      if (same_line(e, L)) {
        already = true;
        break;
      }
    }
    if (already) continue;
    used.push_back(L);

    std::vector<X_monotone_curve_2> pieces;
    for (auto fit = local.faces_begin(); fit != local.faces_end(); ++fit) {
      if (fit->is_fictitious() || fit->data() != FACE_RED) continue;
      clip_line_to_face(fit, L, pieces);
    }
    if (pieces.empty()) continue;

    for (const X_monotone_curve_2& cv : pieces) {
      if (cv.is_segment()) {
        eval_hit(cv.segment().source());
        eval_hit(cv.segment().target());
      } else if (cv.is_ray()) {
        eval_hit(cv.ray().source());
      }
    }
    if (infinite) return;

    CGAL::insert(local, pieces.begin(), pieces.end());
    ++stats.lines_introduced;
    recolor_local_faces(local, M, N, best, consider_F, infinite);
    if (infinite) return;
    merge_green_green(local, stats);
  }
}

// Section 5: d_match = max F. Sample A0, prune green bounded faces with the
// KN variation bound, then incrementally insert A \ A0 on the faces that
// stay red. Returns infinite_F() if some slice with s > 0 has unmatched
// infinite bars.
inline FT maximize_F(const Arrangement_2& arr,
                     const Module& M,
                     const Module& N,
                     bool refine = true,
                     ArrangementStats* stats = nullptr) {
  ArrangementStats local_stats;
  ArrangementStats& st = stats ? *stats : local_stats;
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

  // Vertices of A0 with s ≥ 0 (s = 0 is only used when F is finite).
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    if (vit->is_at_open_boundary()) continue;
    const Point_2 p = vit->point();
    if (p.x() < FT(0)) continue;
    consider_F(p.x(), evaluate_F(M, N, p.x(), p.y()));
    if (infinite) return infinite_F();
  }

  // Midpoints of finite A0 edges. Cheap extra samples on the refine path;
  // the maximizer of F on an A0 edge can sit between vertices (ell vs span
  // at s = 1, t = 1/2). `--no-refine` skips them, so it can still miss.
  if (refine) {
    for (auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
      if (eit->source()->is_at_open_boundary()) continue;
      if (eit->target()->is_at_open_boundary()) continue;
      const Point_2 m = CGAL::midpoint(eit->source()->point(),
                                      eit->target()->point());
      if (m.x() < FT(0)) continue;
      consider_F(m.x(), evaluate_F(M, N, m.x(), m.y()));
      if (infinite) return infinite_F();
    }
  }

  CGAL::Arr_naive_point_location<Arrangement_2> pl(arr);
  std::vector<LabeledGrade> grades;
  collect_labeled_grades(M, N, grades);

  struct FaceWork {
    Arrangement_2::Face_const_handle f;
    Point_2 interior;
    std::vector<Witness> seeds;
    std::vector<Point_2> verts;
    FT F_center;
    FT ub;
    bool unbounded;
    bool has_ub;
  };
  std::vector<FaceWork> works;

  for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    if (fit->is_fictitious()) continue;
    Point_2 interior;
    if (!face_interior_point(pl, fit, interior)) continue;
    if (interior.x() <= FT(0)) continue;

    FaceWork work;
    work.f = fit;
    work.interior = interior;
    work.unbounded = fit->is_unbounded();
    face_finite_vertices(fit, work.verts);

    Witness w0;
    work.F_center = evaluate_F(M, N, interior.x(), interior.y(), &w0);
    consider_F(interior.x(), work.F_center);
    if (infinite) return infinite_F();
    if (w0.valid) work.seeds.push_back(w0);

    // Extra interior samples on bounded cells, refine path only. The
    // centroid of a convex A0 face is inside; one more convex combination
    // of vertices is a cheap extra (λ, μ)-style probe.
    if (refine && !work.unbounded && work.verts.size() >= 3) {
      Point_2 c;
      if (face_centroid(work.verts, c) && point_in_convex_face(fit, c)) {
        Witness wc;
        consider_F(c.x(), evaluate_F(M, N, c.x(), c.y(), &wc));
        if (infinite) return infinite_F();
        if (wc.valid) work.seeds.push_back(wc);
      }
      const Point_2& a = work.verts.front();
      const Point_2& b = work.verts[work.verts.size() / 2];
      const Point_2 p((a.x() + a.x() + b.x()) / FT(3),
                      (a.y() + a.y() + b.y()) / FT(3));
      if (p.x() > FT(0) && point_in_convex_face(fit, p)) {
        Witness wp;
        consider_F(p.x(), evaluate_F(M, N, p.x(), p.y(), &wp));
        if (infinite) return infinite_F();
        if (wp.valid) work.seeds.push_back(wp);
      }
    }

    // s → +∞ along each right-unbounded outer ray of this A0 face.
    // Do not use the paper's c|r + q_x| for the interior witness: on A0
    // that pair generally stops realizing F before infinity, and the
    // closed form then overestimates (often a large integer |q_x − v_x|).
    // The limiting bottleneck is one extra persistence computation per
    // ray, the same cost class as --no-refine.
    if (work.unbounded && interior.x() >= FT(1)) {
      const std::vector<FT> slopes = right_unbounded_slopes(fit);
      for (const FT& r : slopes) {
        consider_F(interior.x(),
                   evaluate_F_s_infinity(M, N, r, interior.x(), interior.y()));
        if (infinite) return infinite_F();
      }
    }

    work.has_ub = false;
    if (!work.unbounded && !is_infinite_F(work.F_center))
      work.has_ub = kn_face_upper_bound(M, N, work.F_center, work.verts, work.ub);
    works.push_back(work);
  }

  if (!refine) return any ? best : FT(0);

  for (FaceWork& work : works) {
    if (infinite) return infinite_F();
    // Green: U(f) ≤ δ, so f cannot beat the samples already in hand.
    // Unbounded A0 faces have no finite KN bound and stay red.
    if (work.has_ub && !(best < work.ub)) continue;
    if (work.seeds.empty()) continue;
    refine_red_a0_face(work.f, work.interior, work.seeds, M, N, grades, pl,
                       work.has_ub, work.ub, consider_F, best, infinite, st);
  }

  return any ? best : FT(0);
}

inline FT maximize_F(const Arrangement_2& arr,
                     const Presentation& M,
                     const Presentation& N,
                     bool refine = true,
                     ArrangementStats* stats = nullptr) {
  return maximize_F(arr, as_module(M), as_module(N), refine, stats);
}

}  // namespace matching
