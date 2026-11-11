#pragma once

// Arrangement with a color per face, for the incremental red/green
// insertion of the A0 segments of one strip (flat or steep):
//   red     F could still exceed Delta somewhere in the face
//   yellow  split from a red face, not classified yet
//   green   F is at most Delta everywhere in the face; never looked at again
//
// The arrangement starts as the quadrilateral that contains every dual
// segment of the strip (see build_initial_quadrilateral).

#include <CGAL/Arr_extended_dcel.h>

#include <matching_distance.h>

namespace matching {

enum class Color { red, yellow, green };

typedef CGAL::Arr_face_extended_dcel<Traits_2, Color> Colored_dcel;
typedef CGAL::Arrangement_2<Traits_2, Colored_dcel> Colored_arrangement;

// Extreme endpoints of the dual segments t = v_y - s v_x, 0 <= s <= 1:
// at s = 0 the segment of v is at t = v_y, at s = 1 at t = v_y - v_x.
struct StripBounds {
  FT t0_min, t0_max;  // on the wall s = 0
  FT t1_min, t1_max;  // on the wall s = 1
};

inline StripBounds strip_bounds(const std::set<Grade2>& grades) {
  assert(!grades.empty());
  const Grade2& g = *grades.begin();
  StripBounds b{g.y, g.y, g.y - g.x, g.y - g.x};
  for (const Grade2& v : grades) {
    b.t0_min = std::min(b.t0_min, v.y);
    b.t0_max = std::max(b.t0_max, v.y);
    b.t1_min = std::min(b.t1_min, v.y - v.x);
    b.t1_max = std::max(b.t1_max, v.y - v.x);
  }
  return b;
}

// Inserts the quadrilateral with corners (0, t0_min), (1, t1_min),
// (1, t1_max), (0, t0_max) into an empty arrangement and colors its
// interior red and the outer face green.
//
// Every dual segment goes from (0, a) to (1, b) with a in [t0_min, t0_max]
// and b in [t1_min, t1_max], so it lies in this quadrilateral. The outer
// face therefore never contains a vertex of A0 and is never split. The
// four corners are endpoints of dual segments, so they are vertices of A0
// themselves.
//
// Degenerate cases: if all endpoints on one wall coincide, that wall is a
// single point and the region is a triangle. If they coincide on both walls,
// all grades have the same dual segment, the region is that segment, and
// there is no red face. Returns the number of red faces (0 or 1).
inline std::size_t build_initial_quadrilateral(Colored_arrangement& arr,
    const Module& M, const Module& N) {
  assert(arr.is_empty());
  std::set<Grade2> grades;
  collect_grades(M, N, grades);
  if (grades.empty()) return 0;
  const StripBounds b = strip_bounds(grades);

  const Point_2 bottom_left(FT(0), b.t0_min), top_left(FT(0), b.t0_max);
  const Point_2 bottom_right(FT(1), b.t1_min), top_right(FT(1), b.t1_max);

  std::vector<X_monotone_curve_2> sides;
  sides.push_back(X_monotone_curve_2(bottom_left, bottom_right));
  if (b.t0_min < b.t0_max || b.t1_min < b.t1_max)  // otherwise top = bottom
    sides.push_back(X_monotone_curve_2(top_left, top_right));
  if (b.t0_min < b.t0_max) sides.push_back(X_monotone_curve_2(bottom_left, top_left));
  if (b.t1_min < b.t1_max) sides.push_back(X_monotone_curve_2(bottom_right, top_right));
  CGAL::insert(arr, sides.begin(), sides.end());

  std::size_t red = 0;
  for (auto f = arr.faces_begin(); f != arr.faces_end(); ++f) {
    if (f->is_unbounded()) {
      f->set_data(Color::green);
    } else {
      f->set_data(Color::red);
      ++red;
    }
  }
  return red;
}

}  // namespace matching
