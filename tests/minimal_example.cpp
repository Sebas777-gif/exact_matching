// Small hand-checked cases for matching::slice_barcode and matching::bottleneck_dist.
// Prints each barcode and exits with 1 if any case is wrong.

#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <colored_arrangement.h>
#include <matching_distance.h>

using namespace matching;

namespace {

Presentation make_presentation(std::vector<Grade2> gens,
    std::vector<std::pair<Grade2, std::vector<int> > > rels) {
  Presentation P;
  P.generators = std::move(gens);
  for (auto& r : rels) {
    P.relations.push_back(r.first);
    P.columns.push_back(r.second);
  }
  return P;
}

// Infinite bars sort after finite ones with the same birth.
bool bar_less(const Bar& a, const Bar& b) {
  if (a.birth != b.birth) return a.birth < b.birth;
  if (a.infinite != b.infinite) return b.infinite;
  return !a.infinite && a.death < b.death;
}

bool bar_equal(const Bar& a, const Bar& b) {
  return a.birth == b.birth && a.infinite == b.infinite && (a.infinite || a.death == b.death);
}

Bar fin(FT birth, FT death) { return Bar{birth, death, false}; }
Bar inf(FT birth) { return Bar{birth, FT(0), true}; }

void print(const std::vector<Bar>& bars) {
  for (const Bar& b : bars) {
    std::cout << " [" << CGAL::to_double(b.birth) << ", ";
    if (b.infinite) std::cout << "inf)";
    else std::cout << CGAL::to_double(b.death) << ")";
  }
  std::cout << "\n";
}

int failures = 0;

void check(const char* name, const Presentation& P, FT s, FT t, std::vector<Bar> expected) {
  std::vector<Bar> got = slice_barcode(P, s, t);
  std::sort(got.begin(), got.end(), bar_less);
  std::sort(expected.begin(), expected.end(), bar_less);
  bool ok = got.size() == expected.size()
      && std::equal(got.begin(), got.end(), expected.begin(), bar_equal);
  std::cout << (ok ? "ok    " : "FAIL  ") << name << "\n  got:     ";
  print(got);
  if (!ok) {
    std::cout << "  expected:";
    print(expected);
    failures++;
  }
}

void print(const std::optional<FT>& d) {
  if (d) std::cout << CGAL::to_double(*d) << "\n";
  else std::cout << "inf\n";
}

void check_bottleneck(const char* name, const std::vector<Bar>& b1, const std::vector<Bar>& b2,
                      const std::optional<FT>& expected) {
  std::optional<FT> got = bottleneck_dist(b1, b2);
  std::optional<FT> swapped = bottleneck_dist(b2, b1);
  bool ok = got.has_value() == expected.has_value() && (!got || *got == *expected)
      && swapped.has_value() == got.has_value() && (!got || *swapped == *got);
  std::cout << (ok ? "ok    " : "FAIL  ") << "bottleneck: " << name << "\n  got:      ";
  print(got);
  if (!ok) {
    std::cout << "  swapped:  ";
    print(swapped);
    std::cout << "  expected: ";
    print(expected);
    failures++;
  }
}

// Checks the initial quadrilateral of the flat strip of (M, N): its counts
// (outer face included), that the bounded faces are red and the outer face
// green, and that inserting all A0 segments into it adds no vertices
// besides those of A0.
void check_quadrilateral(const char* name, const Module& M, const Module& N,
                         std::size_t vertices, std::size_t edges, std::size_t faces,
                         std::size_t red) {
  Colored_arrangement arr;
  FT delta(0);
  std::size_t got_red = build_initial_quadrilateral(arr, M, N, delta);
  bool colors_ok = true;
  std::size_t red_faces = 0;
  for (auto f = arr.faces_begin(); f != arr.faces_end(); ++f) {
    if (f->is_unbounded()) {
      colors_ok = colors_ok && f->data() == Color::green;
    } else {
      colors_ok = colors_ok && f->data() == Color::red;
      ++red_faces;
    }
  }

  // The quadrilateral plus the A0 segments must have exactly the vertices of A0.
  std::set<Grade2> grades;
  collect_grades(M, N, grades);
  std::vector<X_monotone_curve_2> segments;
  for (const Grade2& g : grades) segments.push_back(X_monotone_curve_2(dual_segment_of_grade(g)));
  CGAL::insert(arr, segments.begin(), segments.end());
  Arrangement_2 a0;
  build_arrangement_A(a0, M, N);

  Colored_arrangement fresh;
  FT fresh_delta(0);
  build_initial_quadrilateral(fresh, M, N, fresh_delta);
  bool ok = fresh.number_of_vertices() == vertices && fresh.number_of_edges() == edges
      && fresh.number_of_faces() == faces && got_red == red && red_faces == red && colors_ok
      && arr.number_of_vertices() == a0.number_of_vertices();
  std::cout << (ok ? "ok    " : "FAIL  ") << "quadrilateral: " << name << "\n  got:       "
            << fresh.number_of_vertices() << " vertices, " << fresh.number_of_edges()
            << " edges, " << fresh.number_of_faces() << " faces, " << got_red << " red"
            << (colors_ok ? "" : ", wrong colors") << "; with A0 segments "
            << arr.number_of_vertices() << " vertices (A0: " << a0.number_of_vertices() << ")\n";
  if (!ok) failures++;
}

void check(const char* name, bool ok, const std::string& got) {
  std::cout << (ok ? "ok    " : "FAIL  ") << name << "\n  got:       " << got << "\n";
  if (!ok) failures++;
}

std::string str(const std::optional<FT>& v) {
  if (!v) return "inf";
  std::ostringstream out;
  out << CGAL::exact(*v);
  return out.str();
}

}  // namespace

int main() {
  // g0 = (5,5), g1 = (1,1), g2 = (3,3) in every case.
  const std::vector<Grade2> gens = { {5, 5}, {1, 1}, {3, 3} };

  // The Hera minimal example: r0 = g0 at (6,6).
  // On y = x the sort order is (1, 2, 0), which is not its own inverse.
  Presentation P1 = make_presentation(gens, { { {6, 6}, {0} } });
  check("one relation, y = x", P1, FT(1), FT(0),
        { fin(5, 6), inf(1), inf(3) });

  // Same module on y = x/2 + 1. g0, g2 and r0 lie above the slice and
  // push to their y-coordinate; g1 lies below and pushes to 1/2 + 1.
  check("one relation, y = x/2 + 1", P1, FT(1) / 2, FT(1),
        { fin(5, 6), inf(FT(3) / 2), inf(3) });

  // Same module on the steep slice y = 2x - 4, which is x = y/2 + 2: the
  // slice (s, t) = (1/2, 2) of the swapped chart. There every grade pushes
  // to max(x, y/2 + 2): g0 = (5,5) to 5, g1 = (1,1) to 5/2, g2 = (3,3) to
  // 7/2, and r0 = (6,6) to 6.
  Presentation P1_swapped = swapped_module(as_module(P1)).summands[0];
  check("one relation, y = 2x - 4 (swapped chart)", P1_swapped, FT(1) / 2, FT(2),
        { fin(5, 6), inf(FT(5) / 2), inf(FT(7) / 2) });

  // Relations listed out of order: r0 = g1 at (7,7), r1 = g1 + g2 at (4,4).
  // On y = x, r1 comes first and kills the younger g2 at 4. Then g2 = g1,
  // and r0 kills g1 at 7.
  Presentation P2 = make_presentation(gens, { { {7, 7}, {1} }, { {4, 4}, {1, 2} } });
  check("two relations, out of order", P2, FT(1), FT(0),
        { fin(3, 4), fin(1, 7), inf(5) });

  // r0 = g2 at (3,3) kills g2 at its own grade: a bar of length zero,
  // which must not be reported.
  Presentation P3 = make_presentation(gens, { { {3, 3}, {2} } });
  check("zero-length bar", P3, FT(1), FT(0),
        { inf(1), inf(5) });

  // No relations: every generator is an infinite bar.
  Presentation P4 = make_presentation(gens, {});
  check("free module", P4, FT(1), FT(0),
        { inf(1), inf(3), inf(5) });

  // Bottleneck distances. bottleneck_dist returns std::nullopt for infinity.
  check_bottleneck("both empty", {}, {}, FT(0));
  check_bottleneck("one bar to the diagonal", { fin(0, 4) }, {}, FT(2));
  check_bottleneck("shifted birth", { fin(0, 4) }, { fin(1, 4) }, FT(1));
  // Pairing [0,1) with [0,3) costs 2; sending both to the diagonal costs 3/2.
  check_bottleneck("diagonal beats pairing", { fin(0, 1) }, { fin(0, 3) }, FT(3) / 2);
  // [0,10) pairs with [1,10) at cost 1, the short [3,4) goes to the diagonal at 1/2.
  check_bottleneck("pair plus diagonal", { fin(0, 10) }, { fin(1, 10), fin(3, 4) }, FT(1));
  // Pairing [0,4) with [3,7) costs 3; pairing in order costs 1.
  check_bottleneck("crossing pairs", { fin(0, 4), fin(2, 6) }, { fin(3, 7), fin(1, 5) }, FT(1));
  // Infinite bars are paired by birth: |1 - 2| and |3 - 6|.
  check_bottleneck("infinite bars", { inf(3), inf(1) }, { inf(2), inf(6) }, FT(3));
  // The finite part ([0,10) to the diagonal, 5) dominates the infinite part (1).
  check_bottleneck("infinite and finite", { inf(1), fin(0, 10) }, { inf(2) }, FT(5));
  check_bottleneck("different numbers of infinite bars", { inf(1), inf(3) }, { inf(1) }, std::nullopt);

  // Initial quadrilateral of the flat strip.
  // ell vs span: grades (1,0), (1,2), (0,2), (3,2), (0,0). Endpoints at s = 0
  // (t = y) span [0, 2], at s = 1 (t = y - x) span [-1, 2]: a quadrilateral.
  {
    Presentation ell = make_presentation({ {1, 0}, {0, 2} },
        { { {1, 2}, {0} }, { {3, 2}, {1} } });
    Presentation span = make_presentation({ {0, 0} }, { { {3, 2}, {0} } });
    check_quadrilateral("ell vs span", as_module(ell), as_module(span), 4, 4, 2, 1);
  }
  // Intervals: grades (0,0), (1,1), (2,2). At s = 1 all endpoints are
  // t = 0, so the wall s = 1 is a point and the region is a triangle.
  {
    Presentation long_iv = make_presentation({ {0, 0} }, { { {2, 2}, {0} } });
    Presentation short_iv = make_presentation({ {0, 0} }, { { {1, 1}, {0} } });
    check_quadrilateral("intervals (triangle)", as_module(long_iv), as_module(short_iv),
                        3, 3, 2, 1);
  }
  // One grade (0,0) on both sides: the region is the single segment t = 0,
  // with no red face.
  {
    Presentation free_one = make_presentation({ {0, 0} }, {});
    check_quadrilateral("single grade (segment)", as_module(free_one), as_module(free_one),
                        2, 1, 1, 0);
  }

  // The bar-box bound on the worked example of notes/bar_box_bound: two
  // generators at (0,0) and one relation, their sum, at (50,0) for M and at
  // (0,50) for N. On the triangle T with corners (0,0), (1,0), (1,-50) of the
  // flat strip, F = min(50 - w, 25) with w = 50 s + t, so sup_T F = 25 on a
  // plateau. The variation bound gives 50 there; the bar-box bound is exact.
  {
    Presentation pm = make_presentation({ {0, 0}, {0, 0} }, { { {50, 0}, {0, 1} } });
    Presentation pn = make_presentation({ {0, 0}, {0, 0} }, { { {0, 50}, {0, 1} } });
    const BarBoxModule M = bar_box_module(as_module(pm)), N = bar_box_module(as_module(pn));
    const std::vector<Point_2> T = { Point_2(0, 0), Point_2(1, 0), Point_2(1, -50) };
    IntervalArithmetic arith;
    std::vector<Interval> s, t;
    bar_box_corners<IntervalArithmetic>(T, s, t);
    const double U = bar_box_upper(bar_boxes(M, s, t, arith).boxes, bar_boxes(N, s, t, arith).boxes);
    check("bar-box bound: worked example, upper bound 25", U == 25, std::to_string(U));
    // Delta = 25 is a tie: decided exactly.
    const bool green25 = bar_box_green(M, N, T, FT(25));
    const bool green_below = bar_box_green(M, N, T, FT(25) - FT(1) / 1000);
    check("bar-box bound: green at Delta = 25, red just below", green25 && !green_below,
          std::string(green25 ? "green" : "red") + " / " + (green_below ? "green" : "red"));
  }

  // The whole algorithm with the local phase: matching distances.
  {
    IncrementalOptions options;
    options.local_phase = true;
    auto distance = [&](const Presentation& a, const Presentation& b) {
      return incremental_max_F(as_module(a), as_module(b), options);
    };
    // Intervals [(0,0), (2,2)) and [(0,0), (1,1)): d = 1, attained at an A0 vertex.
    IncrementalResult r1 = distance(make_presentation({ {0, 0} }, { { {2, 2}, {0} } }),
                                    make_presentation({ {0, 0} }, { { {1, 1}, {0} } }));
    check("distance: intervals", r1.distance && *r1.distance == 1 && r1.local.over_limit == 0,
          str(r1.distance));
    // The L of two intervals vs one interval: max F on the A0 vertices is 1,
    // but d = 3/2 at (s, t) = (1, 1/2), on an edge of A0.
    IncrementalResult r2 = distance(
        make_presentation({ {1, 0}, {0, 2} }, { { {1, 2}, {0} }, { {3, 2}, {1} } }),
        make_presentation({ {0, 0} }, { { {3, 2}, {0} } }));
    check("distance: ell vs span (beyond the A0 vertices)",
          r2.a0_max && *r2.a0_max == 1 && r2.distance && *r2.distance == FT(3) / 2
              && r2.local.over_limit == 0,
          "A0 " + str(r2.a0_max) + ", distance " + str(r2.distance));
    // The same pair with at most 2 candidates per face: its faces are split.
    // The pieces around the maximum keep 2 or more candidates, so the result
    // is a lower bound, but it must still reach 3/2 and not exceed it.
    {
      IncrementalOptions split = options;
      split.max_candidates = 2;
      split.split_depth = 3;
      IncrementalResult r = incremental_max_F(
          as_module(make_presentation({ {1, 0}, {0, 2} }, { { {1, 2}, {0} }, { {3, 2}, {1} } })),
          as_module(make_presentation({ {0, 0} }, { { {3, 2}, {0} } })), split);
      check("distance: ell vs span with split faces", r.distance && *r.distance == FT(3) / 2
                && r.local.splits > 0,
            str(r.distance) + ", " + std::to_string(r.local.splits) + " splits");
    }
    // One vs two free generators: different numbers of infinite bars.
    IncrementalResult r3 = distance(make_presentation({ {0, 0} }, {}),
                                    make_presentation({ {0, 0}, {1, 1} }, {}));
    check("distance: unmatched infinite bars", !r3.distance, str(r3.distance));
  }

  if (failures) std::cout << failures << " case(s) failed\n";
  else std::cout << "all cases passed\n";
  return failures ? 1 : 0;
}
