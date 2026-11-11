// Small hand-checked cases for matching::slice_barcode and matching::bottleneck_dist.
// Prints each barcode and exits with 1 if any case is wrong.

#include <iostream>
#include <optional>
#include <tuple>
#include <vector>

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

  // Same module on the steep slice y = 2x - 4. g1 and g2 lie above it and
  // push to (y + 4) / 2; g0 and r0 lie below it and push to x.
  check("one relation, y = 2x - 4", P1, FT(2), FT(-4),
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

  if (failures) std::cout << failures << " case(s) failed\n";
  else std::cout << "all cases passed\n";
  return failures ? 1 : 0;
}
