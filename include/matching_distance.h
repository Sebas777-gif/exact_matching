#pragma once

// Coarse arrangement A0 in the (s, t) plane of slices y = s x + t.
// One dual ray per distinct input grade, plus the walls s = 0 and s = 1.

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <cassert>
#include <set>
#include <utility>
#include <vector>

namespace matching {

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef CGAL::Arr_linear_traits_2<Kernel> Traits_2;
typedef Traits_2::Point_2 Point_2;  // (s, t)
typedef Traits_2::Line_2 Line_2;    // s = 0 and s = 1
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

// Finite presentation. columns is the boundary matrix; A0 uses only the grades.
struct Presentation {
  std::vector<Grade2> generators;
  std::vector<Grade2> relations;
  std::vector<std::vector<int> > columns;
};

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
  if (p.y > s * p.x + t) // strictly above the slice
  {
    if (s >= 1)
      return (p.y - t) / s;
    return p.y;
  }
  if (s >= 1)
    return p.x;
  return s * p.x + t;
}

}  // namespace matching
