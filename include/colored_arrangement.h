#pragma once

// Arrangement with a color per face, for the incremental red/green
// insertion of the A0 segments of one strip (flat or steep):
//   red     F could still exceed Delta somewhere in the face
//   yellow  split from a red face, not classified yet
//   green   F is at most Delta everywhere in the face; never looked at again
//
// The arrangement starts as the quadrilateral that contains every dual
// segment of the strip (see build_initial_quadrilateral).
//
// After all segments are inserted, the red faces are the only places where F
// can exceed the maximum over the A0 vertices. The local phase (local_phase,
// at the end of this file) examines them and gives the matching distance.

#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_observer.h>

#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <random>

#include <bar_box_bound.h>
#include <matching_distance.h>

namespace matching {

enum class Color { red, yellow, green };

// Vertex data: an upper bound on F at the vertex, exact where F was
// evaluated; empty until the vertex is first needed. Halfedge data: unused.
// Face data: the color.
typedef CGAL::Arr_extended_dcel<Traits_2, std::optional<FT>, bool, Color> Colored_dcel;
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

// The constant X of the face bound: the largest |v_x| over all grades v of
// M and N. Moving the slice from (s, t) to (s', t') changes every push
// max(v_y, s v_x + t) by at most X |s - s'| + |t - t'|. For the steep strip,
// pass the swapped modules: their x-coordinates are the original y's.
inline FT max_abs_x(const Module& M, const Module& N) {
  std::set<Grade2> grades;
  collect_grades(M, N, grades);
  FT X(0);
  for (const Grade2& v : grades) X = std::max(X, CGAL::abs(v.x));
  return X;
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
    const Module& M, const Module& N, FT& delta) {
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
  for (auto f = arr.faces_begin(); f != arr.faces_end(); ++f)
    if (f->is_unbounded())
      f->set_data(Color::green);
    else {
      f->set_data(Color::red);
      ++red;
    }
  for (auto v = arr.vertices_begin(); v != arr.vertices_end(); ++v) {
    std::optional<FT> Fv = F(M, N, v->point().x(), v->point().y());
    v->set_data(Fv);
    delta = std::max(delta, *Fv);
  }
  return red;
}

// Which upper bound on F decides whether a face turns green.
//   lipschitz  F(c) + 2 max_i (X |s_i - s_c| + |t_i - t_c|), with one X for
//              all grades (face_upper_bound_lipschitz)
//   variation  F(c) + V_M(c) + V_N(c), with the exact push change of every
//              grade between the corners (face_upper_bound_variation)
//   barbox     the bar-box bound (bar_box_bound.h), which needs no values of F
//   both       variation first (cheap), then barbox
enum class FaceBound { lipschitz, variation, barbox, both };

// The distinct grades of one module (generators and relations of all
// summands).
inline std::vector<Grade2> module_grades(const Module& M) {
  std::set<Grade2> grades;
  for (const Presentation& P : M.summands) {
    grades.insert(P.generators.begin(), P.generators.end());
    grades.insert(P.relations.begin(), P.relations.end());
  }
  return std::vector<Grade2>(grades.begin(), grades.end());
}

// Grades in interval arithmetic for the variation bound (see
// face_upper_bound_variation). Interval is defined in bar_box_bound.h.
// A grade with interval coordinates, converted once per strip.
struct IntervalGrade {
  Interval x, y;
};

inline std::vector<IntervalGrade> interval_grades(const std::vector<Grade2>& grades) {
  std::vector<IntervalGrade> out;
  out.reserve(grades.size());
  for (const Grade2& v : grades)
    out.push_back(IntervalGrade{Interval(CGAL::to_interval(v.x)), Interval(CGAL::to_interval(v.y))});
  return out;
}

// A slice (s, t) where F was evaluated before the insertion.
struct Sample {
  Point_2 point;
  FT F;
};

// The smallest upper bound on F at q that the samples give:
//   F(q) <= F(p) + 2 (X |s_q - s_p| + |t_q - t_p|)   for every sample p,
// the same bound as in face_upper_bound. Empty if there are no samples.
inline std::optional<FT> sample_bound(const std::vector<Sample>& samples, const FT& max_x,
                                      const Point_2& q) {
  std::optional<FT> best;
  for (const Sample& p : samples) {
    const FT bound = p.F + 2 * (max_x * CGAL::abs(q.x() - p.point.x())
                                + CGAL::abs(q.y() - p.point.y()));
    if (!best || bound < *best) best = bound;
  }
  return best;
}

// Watches a Colored_arrangement while segments are inserted and reacts to
// every face split. CGAL calls these functions from inside CGAL::insert.
//
// Usage:
//   Coloring_observer obs(arr, M, N, delta);   // attaches itself to arr
//   CGAL::insert(arr, X_monotone_curve_2(segment));
//    now obs.yellow holds the new faces, obs.green_edges the edges to remove
//
// Do not change the arrangement inside a callback (no remove_edge, no
// insert): CGAL is in the middle of an update. Collect what to do and do
// it after CGAL::insert has returned.
class Coloring_observer : public CGAL::Arr_observer<Colored_arrangement> {
 public:
  typedef Colored_arrangement::Face_handle Face_handle;
  typedef Colored_arrangement::Halfedge_handle Halfedge_handle;

  // Faces split from a red face, still to be classified.
  std::vector<Face_handle> yellow;
  // New edges that split a green face. Both sides are green, so they can be
  // removed once the insertion is finished.
  std::vector<Halfedge_handle> green_edges;
  // Number of F evaluations done by the callbacks.
  std::size_t evaluations = 0;
  // Number of vertices that got a sample bound instead of an evaluation.
  std::size_t skipped = 0;

  // With samples, a new vertex whose sample bound is at most delta is not
  // evaluated: F <= delta there, so it cannot raise delta. The bound is
  // stored instead, which is all face_upper_bound needs. samples must
  // outlive the observer.
  // bound selects the face bound that insert_segment uses.
  Coloring_observer(Colored_arrangement& arr, const Module& M, const Module& N, FT& delta,
                    const std::vector<Sample>* samples = nullptr,
                    FaceBound bound = FaceBound::both)
      : CGAL::Arr_observer<Colored_arrangement>(arr), M_(M), N_(N), delta_(delta),
        max_x_(max_abs_x(M, N)), samples_(samples), bound_(bound),
        grades_M_(interval_grades(module_grades(M))),
        grades_N_(interval_grades(module_grades(N))),
        bar_box_M_(bar_box_module(M)), bar_box_N_(bar_box_module(N)) {}

  // X of the Lipschitz bound for this strip (see max_abs_x).
  const FT& max_x() const { return max_x_; }
  FaceBound bound() const { return bound_; }
  // The distinct grades of each module, as intervals, for the variation bound.
  const std::vector<IntervalGrade>& interval_grades_M() const { return grades_M_; }
  const std::vector<IntervalGrade>& interval_grades_N() const { return grades_N_; }
  // The summands of each module, for the bar-box bound.
  const BarBoxModule& bar_box_M() const { return bar_box_M_; }
  const BarBoxModule& bar_box_N() const { return bar_box_N_; }

  // Called before face f is split by the new edge e. f still has its color.
  // e->source()->point() and e->target()->point() are the two points (s, t)
  // where the new segment crosses the boundary of f.
  void before_split_face(Face_handle f, Halfedge_handle e) override {
    assert(f->data() != Color::yellow); // the queue is cleared after every insertion
    if (f->data() == Color::red) {
      set_value(e->source());
      set_value(e->target());
      std::optional<FT> Fp = e->source()->data();
      std::optional<FT> Fq = e->target()->data();
      assert(Fp && Fq); // F is finite everywhere or nowhere; checked before inserting
      // A stored sample bound is at most delta_, so only evaluated values can
      // raise it.
      delta_ = std::max({*Fp, *Fq, delta_});
    }
    if (f->data() == Color::green)
      green_edges.push_back(e);
  }

  // Called after the split. f is now one half, new_f the other. new_f's
  // color is not copied from f, so set it here.
  void after_split_face(Face_handle f, Face_handle new_f, bool is_hole) override {
    // Every segment goes from the wall s = 0 to the wall s = 1, so a split
    // never creates a hole.
    assert(!is_hole);
    if (f->data() == Color::red) {
      f->set_data(Color::yellow);
      new_f->set_data(Color::yellow);
      yellow.push_back(f);
      yellow.push_back(new_f);
    }
    if (f->data() == Color::green) {
      new_f->set_data(Color::green);
    }
  }

 private:
  // Gives v its value, unless it has one: the sample bound if that is at
  // most delta_ (covered by a sample), otherwise F itself.
  void set_value(Colored_arrangement::Vertex_handle v) {
    if (v->data()) return;
    if (samples_) {
      const std::optional<FT> bound = sample_bound(*samples_, max_x_, v->point());
      if (bound && *bound <= delta_) {
        v->set_data(bound);
        ++skipped;
        return;
      }
    }
    v->set_data(F(M_, N_, v->point().x(), v->point().y()));
    ++evaluations;
  }

  const Module& M_;
  const Module& N_;
  FT& delta_;  // the largest F found so far, shared by both strips
  FT max_x_;   // the largest |v_x| over all grades v of M_ and N_, fixed per strip
  const std::vector<Sample>* samples_;  // samples of this strip, or nullptr
  FaceBound bound_;
  std::vector<IntervalGrade> grades_M_, grades_N_;
  BarBoxModule bar_box_M_, bar_box_N_;
};

// The corners of a bounded face, counterclockwise.
template <class FaceHandle>
std::vector<Point_2> face_corners(FaceHandle f) {
  std::vector<Point_2> corners;
  auto he = f->outer_ccb();
  auto first = he;
  do {
    corners.push_back(he->target()->point());
  } while (++he != first);
  return corners;
}

inline FT face_upper_bound_lipschitz(Colored_arrangement::Face_handle f, const FT& max_x) {
  // All corners of a bounded face f, in counterclockwise order.
  std::vector<Colored_arrangement::Vertex_handle> corners;
  auto he = f->outer_ccb();     // circulator over the halfedges around f
  auto first = he;
  do {
    corners.push_back(he->target());
  } while (++he != first);
  // The bound, for one corner c = (s_c, t_c) of f whose F is known:
  //
  //   Moving the slice from c to any p = (s_p, t_p) moves every push
  //   max(v_y, s v_x + t) by at most  X |s_p - s_c| + |t_p - t_c|.
  //   Both modules move, so F moves by at most twice that:
  //
  //     F(p) <= F(c) + 2 (X |s_p - s_c| + |t_p - t_c|)     for every p.
  //
  //   The right-hand side is convex in p, so on the convex face f it is
  //   largest at one of the corners (s_i, t_i). Hence for every p in f:
  //
  //     F(p) <= F(c) + 2 max_i (X |s_i - s_c| + |t_i - t_c|).
  //
  // Every corner can play c; the smallest of these bounds is returned.
  // (min_bound must not start at delta: that would cap every bound at delta
  // and make every face green.)
  std::optional<FT> min_bound;
  for (Colored_arrangement::Vertex_handle c : corners) {
    assert(c->data());  // every corner of a red or yellow face has F stored
    const FT& s_c = c->point().x();
    const FT& t_c = c->point().y();

    // farthest = max_i (X |s_i - s_c| + |t_i - t_c|): the distance from c to
    // the farthest corner of f, measured so that it bounds the change of F.
    FT farthest(0);
    for (Colored_arrangement::Vertex_handle i : corners) {
      const FT dist = max_x * CGAL::abs(i->point().x() - s_c) + CGAL::abs(i->point().y() - t_c);
      farthest = std::max(farthest, dist);
    }

    const FT bound = *c->data() + 2 * farthest;
    if (!min_bound || bound < *min_bound) min_bound = bound;
  }
  return *min_bound;
}

// Upper bound on F over the face f, from the push changes of the grades
// between the corners of f.
//
// For a corner c of f whose F is known (exact or an upper bound) and any
// slice p in f, let d_v(p) = push_v(p) - push_v(c). By stability, applied to
// each module separately,
//
//   F(p) <= F(c) + max_{v in M} |d_v(p)| + max_{v in N} |d_v(p)|.
//
// push_v = max(v_y, s v_x + t) takes its largest and its smallest value over
// the convex face f at corners: the largest because it is convex; the
// smallest because it is v_y (its global minimum) on the half-plane
// s v_x + t <= v_y, which contains a corner of f if it meets f at all, and
// otherwise push_v is linear on f. So with corners c_i,
//
//   V_M(c) = max_{v in M} max(max_i push_v(c_i) - push_v(c), push_v(c) - min_i push_v(c_i))
//
// bounds the first maximum on all of f, likewise V_N(c), and
//
//   sup_f F <= min over corners c of F(c) + V_M(c) + V_N(c).
//
// A grade that lies above every slice in f (push = v_y throughout)
// contributes nothing, which is where this beats the Lipschitz bound.
//
// The bound only has to be an upper bound, not exact, so it is computed in
// interval arithmetic (doubles rounded outward) instead of exact FT, which
// is much faster. The result is the upper end: a valid upper bound on F.
inline double face_upper_bound_variation(Colored_arrangement::Face_handle f,
                                         const std::vector<IntervalGrade>& grades_M,
                                         const std::vector<IntervalGrade>& grades_N) {
  std::vector<Colored_arrangement::Vertex_handle> corners;
  auto he = f->outer_ccb();
  auto first = he;
  do {
    corners.push_back(he->target());
  } while (++he != first);
  const std::size_t k = corners.size();
  std::vector<Interval> s(k), t(k);
  for (std::size_t i = 0; i < k; i++) {
    s[i] = Interval(CGAL::to_interval(corners[i]->point().x()));
    t[i] = Interval(CGAL::to_interval(corners[i]->point().y()));
  }

  // For one module: the push of every grade at every corner, and per grade
  // an upper bound on its largest and a lower bound on its smallest push
  // over the corners. max is monotone, so the ends of max(a, b) are the
  // maxima of the ends.
  struct Pushes {
    std::vector<std::vector<Interval> > at;  // at[i][v]: push of grade v at corner i
    std::vector<double> hi, lo;              // per grade: >= max and <= min over the corners
  };
  auto pushes = [&](const std::vector<IntervalGrade>& grades) {
    Pushes P;
    P.at.assign(k, std::vector<Interval>(grades.size()));
    for (std::size_t i = 0; i < k; i++)
      for (std::size_t v = 0; v < grades.size(); v++) {
        const Interval line = s[i] * grades[v].x + t[i];
        P.at[i][v] = Interval(std::max(grades[v].y.inf(), line.inf()),
                              std::max(grades[v].y.sup(), line.sup()));
      }
    P.hi.assign(grades.size(), -std::numeric_limits<double>::infinity());
    P.lo.assign(grades.size(), std::numeric_limits<double>::infinity());
    for (std::size_t i = 0; i < k; i++)
      for (std::size_t v = 0; v < grades.size(); v++) {
        P.hi[v] = std::max(P.hi[v], P.at[i][v].sup());
        P.lo[v] = std::min(P.lo[v], P.at[i][v].inf());
      }
    return P;
  };
  // Upper bound on V(c): the largest change of any grade's push between
  // corner c and the rest of f.
  auto variation = [](const Pushes& P, std::size_t c) {
    double V = 0;
    for (std::size_t v = 0; v < P.hi.size(); v++) {
      V = std::max(V, (Interval(P.hi[v]) - P.at[c][v]).sup());
      V = std::max(V, (P.at[c][v] - Interval(P.lo[v])).sup());
    }
    return V;
  };

  const Pushes PM = pushes(grades_M), PN = pushes(grades_N);
  double min_bound = std::numeric_limits<double>::infinity();
  for (std::size_t c = 0; c < k; c++) {
    assert(corners[c]->data());  // every corner of a red or yellow face has a value
    const Interval bound = Interval(CGAL::to_interval(*corners[c]->data()))
        + Interval(variation(PM, c)) + Interval(variation(PN, c));
    min_bound = std::min(min_bound, bound.sup());
  }
  return min_bound;
}

// Whether F <= delta is guaranteed on all of f, with the bound selected in
// the observer. The variation bound is an interval upper end, so it is
// compared with the lower end of delta: rounding can only make a face stay
// red, never wrongly turn green.
//
// The bar-box bound decides exactly where interval arithmetic is unsure
// (bar_box_green), because on plateaus it often equals delta exactly.
inline bool face_is_green(Colored_arrangement::Face_handle f, const Coloring_observer& obs,
                          const FT& delta) {
  const FaceBound bound = obs.bound();
  if (bound == FaceBound::lipschitz)
    return face_upper_bound_lipschitz(f, obs.max_x()) <= delta;
  if (bound == FaceBound::variation || bound == FaceBound::both) {
    if (face_upper_bound_variation(f, obs.interval_grades_M(), obs.interval_grades_N())
        <= CGAL::to_interval(delta).first)
      return true;
    if (bound == FaceBound::variation) return false;
  }
  return bar_box_green(obs.bar_box_M(), obs.bar_box_N(), face_corners(f), delta);
}

// Inserts one segment, then restores the invariant: no yellow faces, and
// no edge with green on both sides.
inline void insert_segment(Colored_arrangement& arr, Coloring_observer& obs,
                           const Segment_2& seg, const FT& delta) {
  CGAL::insert(arr, X_monotone_curve_2(seg));   // the callbacks fill obs.yellow and obs.green_edges

  // Step 1: classify every yellow face as red or green (face_upper_bound vs delta).
  std::vector<Colored_arrangement::Face_handle> new_green_faces;
  for (Colored_arrangement::Face_handle f : obs.yellow)
    if (face_is_green(f, obs, delta)) {
      f->set_data(Color::green);
      new_green_faces.push_back(f);
    }
    else
      f->set_data(Color::red);
  // Step 2: collect the edges with green on both sides, once per edge. Edges
  //         of the outer boundary stay: the later segments end on them.
  auto removable = [](Colored_arrangement::Halfedge_handle he) {
    return he->face()->data() == Color::green && he->twin()->face()->data() == Color::green
        && !he->face()->is_unbounded() && !he->twin()->face()->is_unbounded();
  };
  std::set<Colored_arrangement::Halfedge_handle> to_remove;  // one halfedge per edge
  auto collect = [&](Colored_arrangement::Halfedge_handle he) {
    if (removable(he)) to_remove.insert(&*he < &*he->twin() ? he : he->twin());
  };
  for (Colored_arrangement::Face_handle f : new_green_faces) {
    auto he = f->outer_ccb();
    auto first = he;
    do collect(he); while (++he != first);
  }
  for (Colored_arrangement::Halfedge_handle he : obs.green_edges) collect(he);

  // Step 3: remove them, then clear the observer's lists.
  obs.yellow.clear();
  obs.green_edges.clear();
  for (Colored_arrangement::Halfedge_handle he : to_remove) arr.remove_edge(he);
}

// --- Randomized insertion --------------------------------------------------

// Counts of the incremental algorithm for one strip.
struct StripStats {
  std::size_t segments = 0;     // dual segments inserted
  std::size_t samples = 0;      // F evaluations of the sampling pass
  std::size_t evaluations = 0;  // F evaluations during the insertion, quadrilateral corners included
  std::size_t skipped = 0;      // vertices covered by a sample, not evaluated
  // The final arrangement, outer face included.
  std::size_t vertices = 0, edges = 0, faces = 0, red_faces = 0;
};

// A random vertex of A0 in the strip 0 <= s <= 1: where the dual segments of
// two random grades a and b cross, if they do (a_y - s a_x = b_y - s b_x),
// and otherwise an endpoint of a's segment on a random wall.
//
// Sampling only A0 vertices keeps the result comparable: Delta never
// exceeds the maximum of F over the A0 vertices, so the incremental result
// must equal max_F_on_vertices exactly. (Any slice gives a valid lower
// bound on the matching distance; for the A0 experiment that would break
// the comparison.)
inline Point_2 random_a0_vertex(const std::vector<Grade2>& grades, std::mt19937& rng) {
  std::uniform_int_distribution<std::size_t> pick(0, grades.size() - 1);
  const Grade2& a = grades[pick(rng)];
  const Grade2& b = grades[pick(rng)];
  if (a.x != b.x) {
    const FT s = (b.y - a.y) / (b.x - a.x);
    if (s >= 0 && s <= 1) return Point_2(s, a.y - s * a.x);
  }
  if (rng() % 2) return Point_2(FT(0), a.y);
  return Point_2(FT(1), a.y - a.x);
}

// Evaluates F at k random A0 vertices of the strip of (M, N), raises delta,
// and appends them to samples. Returns false, after the first sample, if F is infinite there:
// F is infinite on every slice or on none. A sampled vertex is evaluated
// again if the insertion reaches it later.
inline bool sample_delta(const Module& M, const Module& N, std::size_t k, std::mt19937& rng,
                         FT& delta, std::vector<Sample>& samples, StripStats& stats) {
  std::set<Grade2> grade_set;
  collect_grades(M, N, grade_set);
  if (grade_set.empty()) return true;
  const std::vector<Grade2> grades(grade_set.begin(), grade_set.end());
  for (std::size_t i = 0; i < k; i++) {
    const Point_2 p = random_a0_vertex(grades, rng);
    const std::optional<FT> f = F(M, N, p.x(), p.y());
    ++stats.samples;
    if (!f) return false;
    delta = std::max(delta, *f);
    samples.push_back(Sample{p, *f});
  }
  return true;
}

// Inserts the dual segments of the strip of (M, N) in random order into the
// initial quadrilateral, with red/green coloring, and raises delta. With
// samples, vertices covered by a sample are not evaluated (see
// Coloring_observer). arr must be empty; afterwards it holds the colored
// arrangement, whose red faces the local phase examines.
inline void insert_strip(Colored_arrangement& arr, const Module& M, const Module& N,
                         std::mt19937& rng, FT& delta, StripStats& stats,
                         const std::vector<Sample>* samples = nullptr,
                         FaceBound bound = FaceBound::both) {
  build_initial_quadrilateral(arr, M, N, delta);
  stats.evaluations = arr.number_of_vertices();  // the corners

  // Detached again when insert_strip returns; arr stays.
  Coloring_observer obs(arr, M, N, delta, samples, bound);

  std::set<Grade2> grades;
  collect_grades(M, N, grades);
  std::vector<Segment_2> segments;
  for (const Grade2& g : grades) segments.push_back(dual_segment_of_grade(g));
  std::shuffle(segments.begin(), segments.end(), rng);
  for (const Segment_2& seg : segments) insert_segment(arr, obs, seg, delta);

  stats.segments = segments.size();
  stats.evaluations += obs.evaluations;
  stats.skipped = obs.skipped;
  stats.vertices = arr.number_of_vertices();
  stats.edges = arr.number_of_edges();
  stats.faces = arr.number_of_faces();
  for (auto f = arr.faces_begin(); f != arr.faces_end(); ++f)
    if (f->data() == Color::red) ++stats.red_faces;
}

// --- Local phase: beyond the A0 vertices -------------------------------------
//
// After all A0 segments are inserted, every red face lies inside one face of
// A0, so every push is linear on it: v_y if v lies on or above the slices of
// the face, s v_x + t otherwise. At every slice, F then equals one of finitely
// many linear functions (KLO): the difference of two births or of two deaths
// of the two modules, or half a bar length within one module.
//
// Let Delta be the largest F found so far and U an upper bound on F over the
// red face f (the bar-box bound). Only functions that can realize a value in
// (Delta, U] on f matter, and the bar-box classes restrict which grades can
// form them ("candidates"). F is continuous and piecewise linear on f, with
// pieces given by the candidates where F > Delta. So if the maximum of F on
// f exceeds Delta, it is attained
//   - at a corner of f (already evaluated),
//   - where two candidates are equal on the boundary of f, or
//   - where three candidates are equal inside f.
// The local phase evaluates F exactly at those points, best candidate value
// first, and raises Delta. No line of the refined arrangement A is inserted.
//
// This argument is in notes/bar_box_bound; it still needs independent
// checking. A face with more than max_candidates candidates is skipped, and
// the result is then only a lower bound (LocalStats::over_limit > 0).

struct LocalStats {
  std::size_t red_faces = 0;    // red faces that the bar-box bound could not turn green
  std::size_t candidates = 0;   // candidate functions, summed over those faces
  std::size_t largest = 0;      // most candidates in one face
  std::size_t lines = 0;        // pairs of candidates whose equality line matters in a face
  std::size_t points = 0;       // points where F was a candidate for evaluation
  std::size_t evaluations = 0;  // F evaluations
  std::size_t certified = 0;    // points shown F <= Delta by the bar-box test, not evaluated
  std::size_t over_limit = 0;   // faces or pieces skipped: too many candidates at max depth
  std::size_t splits = 0;       // faces or pieces cut in two (--split-depth)
  // Level-line certification (--level-lines):
  std::size_t level_tries = 0;      // polygons it was tried on
  std::size_t level_lines = 0;      // lines e = Delta crossing them, summed
  std::size_t level_cells = 0;      // cells of those line arrangements
  std::size_t level_checks = 0;     // cells with a candidate above Delta: F checked at a point
  std::size_t level_certified = 0;  // polygons shown F <= Delta this way
  std::size_t level_raised = 0;     // tries that found F > Delta (then Delta grows)
  std::size_t level_skipped = 0;    // polygons with more level lines than the limit
};

// A linear function a s + b t + c of the slice (s, t).
struct Linear {
  FT a, b, c;
  FT at(const Point_2& p) const { return a * p.x() + b * p.y() + c; }
};

inline bool operator<(const Linear& x, const Linear& y) {
  if (x.a != y.a) return x.a < y.a;
  if (x.b != y.b) return x.b < y.b;
  return x.c < y.c;
}

inline Linear operator-(const Linear& x, const Linear& y) { return {x.a - y.a, x.b - y.b, x.c - y.c}; }

// The same in interval arithmetic, for cheap tests before exact ones.
struct LinearI {
  Interval a, b, c;
  explicit LinearI(const Linear& l)
      : a(CGAL::to_interval(l.a)), b(CGAL::to_interval(l.b)), c(CGAL::to_interval(l.c)) {}
  LinearI(Interval a_, Interval b_, Interval c_) : a(a_), b(b_), c(c_) {}
  Interval at(const Interval& s, const Interval& t) const { return a * s + b * t + c; }
};

inline LinearI operator-(const LinearI& x, const LinearI& y) {
  return LinearI(x.a - y.a, x.b - y.b, x.c - y.c);
}

// The push of v on a face of A0, as a linear function. mid is a point inside
// the face; no dual line crosses the inside of an A0 face.
inline Linear push_on_face(const Grade2& v, const Point_2& mid) {
  if (v.y >= mid.x() * v.x + mid.y()) return {FT(0), FT(0), v.y};
  return {v.x, FT(1), FT(0)};
}

struct PointLess {
  bool operator()(const Point_2& p, const Point_2& q) const {
    return CGAL::compare_xy(p, q) == CGAL::SMALLER;
  }
};

namespace detail {

// The candidate functions of a red face with corners `corners`: every
// function that the bar-box classes allow and that can take a value in
// (delta, U] somewhere on the face.
inline std::vector<Linear> local_candidates(const BarBoxModule& M, const BarBoxModule& N,
                                            const BarBoxes<Interval>& X, const BarBoxes<Interval>& Y,
                                            const std::vector<Point_2>& corners,
                                            const Point_2& mid, const FT& delta, double U) {
  IntervalArithmetic arith;  // only for "may exceed delta" tests; uncertainty keeps a function
  const Interval delta_i = IntervalArithmetic::from(delta);
  std::vector<Interval> s, t;
  bar_box_corners<IntervalArithmetic>(corners, s, t);
  std::map<Linear, char> found;

  // Keeps e if it can take a value in (delta, cap] on the face, cap <= U.
  auto add = [&](const Linear& e, double cap) {
    const LinearI ei(e);
    double hi = -std::numeric_limits<double>::infinity(), lo = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < s.size(); i++) {
      const Interval v = ei.at(s[i], t[i]);
      hi = std::max(hi, v.sup());
      lo = std::min(lo, v.inf());
    }
    // On the face, e is linear, so its extremes are at corners.
    if (hi > delta_i.inf() && lo <= std::min(U, cap)) found.emplace(e, 0);
  };
  auto unique_boxes = [](const BarBoxes<Interval>& B) {
    std::vector<const BarBox<Interval>*> out;
    std::set<std::tuple<std::size_t, int, int> > seen;
    for (const auto& b : B.boxes)
      if (seen.insert(std::make_tuple(b.summand, b.gen_class, b.rel_class)).second) out.push_back(&b);
    return out;
  };
  const auto UX = unique_boxes(X), UY = unique_boxes(Y);

  // Half bar lengths within one module.
  auto half_lengths = [&](const BarBoxModule& mod, const BarBoxes<Interval>& B,
                          const std::vector<const BarBox<Interval>*>& boxes) {
    for (const BarBox<Interval>* b : boxes) {
      if (b->essential || arith.le(bar_box_diag<IntervalArithmetic>(*b), delta_i)) continue;
      const BarBoxSummand& q = mod[b->summand];
      for (int g : B.gen_members[b->summand][b->gen_class])
        for (int r : B.rel_members[b->summand][b->rel_class]) {
          const Linear d = push_on_face(q.rels[r], mid) - push_on_face(q.gens[g], mid);
          add(Linear{d.a / 2, d.b / 2, d.c / 2}, U);
        }
    }
  };
  half_lengths(M, X, UX);
  half_lengths(N, Y, UY);

  // Differences of births, or of deaths, across the modules (both signs).
  //
  // There is always an optimal matching in which every matched pair costs
  // at most the larger half length of its two bars: otherwise sending both
  // bars to the diagonal is cheaper. So a matched pair that realizes F > delta
  // has a bar of half length > delta, and one of its boxes is forced (diag
  // may exceed delta). Pairs of two short bars are left out, and the value of
  // a pair is at most cap = max(diag(x), diag(y)) (infinite for essential bars).
  auto forced = [&](const BarBox<Interval>* b) {
    return b->essential || !arith.le(bar_box_diag<IntervalArithmetic>(*b), delta_i);
  };
  for (const BarBox<Interval>* x : UX)
    for (const BarBox<Interval>* y : UY) {
      if (x->essential != y->essential) continue;
      if (!forced(x) && !forced(y)) continue;
      if (arith.le(bar_box_cost<IntervalArithmetic>(*x, *y), delta_i)) continue;
      const double cap = x->essential ? std::numeric_limits<double>::infinity()
          : std::max(bar_box_diag<IntervalArithmetic>(*x).sup(),
                     bar_box_diag<IntervalArithmetic>(*y).sup());
      const BarBoxSummand& qx = M[x->summand];
      const BarBoxSummand& qy = N[y->summand];
      for (int g : X.gen_members[x->summand][x->gen_class])
        for (int h : Y.gen_members[y->summand][y->gen_class]) {
          const Linear e = push_on_face(qx.gens[g], mid) - push_on_face(qy.gens[h], mid);
          add(e, cap);
          add(Linear{-e.a, -e.b, -e.c}, cap);
        }
      if (x->essential) continue;
      for (int r : X.rel_members[x->summand][x->rel_class])
        for (int q : Y.rel_members[y->summand][y->rel_class]) {
          const Linear e = push_on_face(qx.rels[r], mid) - push_on_face(qy.rels[q], mid);
          add(e, cap);
          add(Linear{-e.a, -e.b, -e.c}, cap);
        }
    }

  std::vector<Linear> out;
  for (const auto& kv : found) out.push_back(kv.first);
  return out;
}

}  // namespace detail

enum class LocalOutcome { green, done, over_limit };

// --- Certifying a polygon via level lines -------------------------------------
//
// Deciding F <= delta on a polygon P inside a face of A0 needs far fewer
// lines than finding the maximum (Bjerkevik and Kerber decide d_M <= lambda
// the same way). Let E be the candidates of P for delta. Wherever F > delta
// in P, F equals some e in E (the argument of the local phase). Then:
//
//   The boundary of {F > delta} inside P lies on the level lines e = delta,
//   e in E. Take p on that boundary and points q_n -> p with F(q_n) > delta.
//   Infinitely many q_n share one candidate e with F(q_n) = e(q_n), and by
//   continuity of F and e, e(p) = F(p) = delta.
//
// So on each cell of the arrangement of the level lines inside P, either
// F > delta everywhere or nowhere. A point of {F > delta} on a line or on the
// edge of P has a neighbourhood in {F > delta}, which meets an open cell, so
// checking one inside point per cell decides F <= delta on all of P. A cell
// where no candidate exceeds delta needs no check at all: F > delta at a
// point means some candidate does. No lines where the order of two pushes
// changes are needed, as F is continuous across them.
//
// The cells are found by cutting P recursively along the level lines (a
// binary space partition), not by building their arrangement. A piece
// needs no further cuts once no candidate exceeds delta on it, or once the
// bar-box bound shows F <= delta on it. Cells are checked exactly at a point.

enum class LevelOutcome { certified, raised, too_many_lines };

namespace detail {

// A convex polygon as its corners and the lines of its edges: edge q, from
// corner q to corner q+1, lies on lines[q]. Corners are always computed as
// the intersection of two edge lines, so they never become nested
// constructions.
struct LinePolygon {
  std::vector<Point_2> corners;
  std::vector<Linear> lines;
};

inline Point_2 meet(const Linear& l1, const Linear& l2) {
  const FT d = l1.a * l2.b - l1.b * l2.a;
  assert(d != 0);
  return Point_2((l1.b * l2.c - l2.b * l1.c) / d, (l2.a * l1.c - l1.a * l2.c) / d);
}

inline Linear line_through(const Point_2& p, const Point_2& q) {
  // (q - p) x (r - p) = 0 as a s + b t + c = 0, positive to the left.
  const FT a = -(q.y() - p.y()), b = q.x() - p.x();
  return Linear{a, b, -(a * p.x() + b * p.y())};
}

// The part of P where L >= 0. L must cross the inside of P.
inline LinePolygon clip_lines(const LinePolygon& P, const Linear& L) {
  LinePolygon out;
  const std::size_t k = P.corners.size();
  std::vector<CGAL::Sign> sg(k);
  for (std::size_t q = 0; q < k; q++) sg[q] = CGAL::sign(L.at(P.corners[q]));
  for (std::size_t q = 0; q < k; q++) {
    const std::size_t r = (q + 1) % k;
    if (sg[q] == CGAL::POSITIVE) {
      out.corners.push_back(P.corners[q]);
      out.lines.push_back(P.lines[q]);
      if (sg[r] == CGAL::NEGATIVE) {
        out.corners.push_back(meet(P.lines[q], L));
        out.lines.push_back(L);
      }
    } else if (sg[q] == CGAL::ZERO) {
      out.corners.push_back(P.corners[q]);
      out.lines.push_back(sg[r] == CGAL::NEGATIVE ? L : P.lines[q]);
    } else if (sg[r] == CGAL::POSITIVE) {
      out.corners.push_back(meet(P.lines[q], L));
      out.lines.push_back(P.lines[q]);
    }
  }
  return out;
}

}  // namespace detail

inline LevelOutcome level_certify(const std::vector<Point_2>& corners, const std::vector<Linear>& E,
                                  const BarBoxModule& MB, const BarBoxModule& NB, const Module& M,
                                  const Module& N, FT& delta, std::size_t max_lines,
                                  LocalStats& stats, std::set<Point_2, PointLess>& evaluated) {
  const std::size_t k = corners.size(), m = E.size();
  // The level line e - delta = 0 of each candidate, scaled so that equal
  // lines coincide; line_of[i] is the index of candidate i's line (or -1
  // for a constant candidate).
  std::map<Linear, int> line_index;
  std::vector<Linear> level;
  std::vector<int> line_of(m, -1);
  for (std::size_t i = 0; i < m; i++) {
    const Linear L{E[i].a, E[i].b, E[i].c - delta};
    if (L.a == 0 && L.b == 0) continue;
    const FT n = L.a != 0 ? L.a : L.b;
    const Linear key{L.a / n, L.b / n, L.c / n};
    auto it = line_index.emplace(key, static_cast<int>(level.size()));
    if (it.second) level.push_back(key);
    line_of[i] = it.first->second;
  }
  if (level.size() > max_lines) {
    ++stats.level_skipped;
    return LevelOutcome::too_many_lines;
  }
  ++stats.level_tries;
  // Cuts in a random order (a randomized binary space partition).
  std::vector<int> rank(level.size());
  std::iota(rank.begin(), rank.end(), 0);
  std::shuffle(rank.begin(), rank.end(), std::mt19937(static_cast<unsigned>(level.size())));

  std::vector<LinearI> EI;
  for (const Linear& e : E) EI.push_back(LinearI(e));
  const Interval delta_i = IntervalArithmetic::from(delta);

  detail::LinePolygon root;
  root.corners = corners;
  for (std::size_t q = 0; q < k; q++)
    root.lines.push_back(detail::line_through(corners[q], corners[(q + 1) % k]));
  std::vector<std::size_t> all(m);
  std::iota(all.begin(), all.end(), 0);
  std::set<int> lines_seen;

  // Returns false once F > delta was found (and delta raised).
  // A candidate that is >= delta on all of a piece and > delta somewhere
  // stays so on every sub-piece, so instead of passing it down, a flag
  // `above` records that one exists. Only candidates whose line crosses the
  // piece are passed on.
  std::function<bool(const detail::LinePolygon&, const std::vector<std::size_t>&, bool, bool)> visit =
      [&](const detail::LinePolygon& P, const std::vector<std::size_t>& parent, bool above,
          bool try_box) -> bool {
    const std::size_t pk = P.corners.size();
    std::vector<Interval> s(pk), t(pk);
    for (std::size_t q = 0; q < pk; q++) {
      s[q] = CGAL::to_interval(P.corners[q].x());
      t[q] = CGAL::to_interval(P.corners[q].y());
    }
    // The candidates of the parent whose level line crosses the inside of P.
    std::vector<std::size_t> C;
    int cut = -1;
    for (std::size_t i : parent) {
      bool up = false, down = false;
      for (std::size_t q = 0; q < pk; q++) {
        const Interval v = EI[i].at(s[q], t[q]);
        if (v.inf() > delta_i.sup()) up = true;
        else if (v.sup() < delta_i.inf()) down = true;
        else {
          const CGAL::Comparison_result c = CGAL::compare(E[i].at(P.corners[q]), delta);
          up = up || c == CGAL::LARGER;
          down = down || c == CGAL::SMALLER;
        }
      }
      if (!up) continue;  // <= delta on P
      if (!down) {        // >= delta on P, > delta inside
        above = true;
        continue;
      }
      C.push_back(i);
      if (cut < 0 || rank[line_of[i]] < rank[cut]) cut = line_of[i];
    }
    if (C.empty() && !above) return true;  // no candidate above delta: F <= delta on P
    // Interval arithmetic only: at plateaus where F = delta, the exact
    // fallback is uncertain on large pieces and costs more than cutting on.
    if (try_box && bar_box_green(MB, NB, P.corners, delta, false)) return true;
    if (cut < 0) {
      // A cell: some candidate exceeds delta on all of it. Check one point.
      ++stats.level_cells;
      ++stats.level_checks;
      FT ms(0), mt(0);
      for (const Point_2& c : P.corners) {
        ms += c.x();
        mt += c.y();
      }
      const Point_2 p(ms / static_cast<int>(pk), mt / static_cast<int>(pk));
      if (!evaluated.insert(p).second) return true;
      if (bar_box_green(MB, NB, std::vector<Point_2>{p}, delta)) {
        ++stats.certified;
        return true;
      }
      const std::optional<FT> value = F(M, N, p.x(), p.y());
      ++stats.evaluations;
      assert(value);
      if (*value > delta) {
        delta = *value;
        return false;
      }
      return true;
    }
    lines_seen.insert(cut);
    const Linear& L = level[cut];
    return visit(detail::clip_lines(P, L), C, above, true)
        && visit(detail::clip_lines(P, Linear{-L.a, -L.b, -L.c}), C, above, true);
  };
  const bool ok = visit(root, all, false, false);  // the whole polygon is not green
  stats.level_lines += lines_seen.size();
  if (!ok) {
    ++stats.level_raised;
    return LevelOutcome::raised;
  }
  ++stats.level_certified;
  return LevelOutcome::certified;
}

// The local phase on one convex polygon of slices that lies inside a single
// face of A0 (a red face, or a piece of one). Raises delta to the maximum of F
// over the polygon, or returns over_limit if it has more than max_candidates
// candidates (then nothing is evaluated apart from new corners). With
// new_corners, F is first checked at the corners, which the insertion has
// not seen.
inline LocalOutcome local_polygon(const std::vector<Point_2>& corners, const BarBoxModule& MB,
                                  const BarBoxModule& NB, const Module& M, const Module& N,
                                  FT& delta, std::size_t max_candidates, bool new_corners,
                                  LocalStats& stats, std::set<Point_2, PointLess>& evaluated,
                                  std::size_t level_lines = 0, std::size_t level_rounds = 3) {
  if (bar_box_green(MB, NB, corners, delta)) return LocalOutcome::green;  // delta may have grown
  const std::size_t k = corners.size();

  // Evaluates F at p unless already done or certified <= delta.
  auto evaluate = [&](const Point_2& p) {
    if (!evaluated.insert(p).second) return;
    if (bar_box_green(MB, NB, std::vector<Point_2>{p}, delta)) {
      ++stats.certified;
      return;
    }
    const std::optional<FT> value = F(M, N, p.x(), p.y());
    ++stats.evaluations;
    assert(value);  // F is finite (checked before)
    if (*value > delta) delta = *value;
  };
  // Corners of a piece of a split face were not seen by the insertion.
  if (new_corners)
    for (const Point_2& c : corners) evaluate(c);

  IntervalArithmetic arith;
  std::vector<Interval> si, ti;
  bar_box_corners<IntervalArithmetic>(corners, si, ti);
  const BarBoxes<Interval> X = bar_boxes(MB, si, ti, arith), Y = bar_boxes(NB, si, ti, arith);
  const double U = bar_box_upper(X.boxes, Y.boxes);
  const bool bounded = U < std::numeric_limits<double>::infinity();
  const FT U_exact = bounded ? FT(U) : FT(0);

  FT ms(0), mt(0);
  for (const Point_2& c : corners) {
    ms += c.x();
    mt += c.y();
  }
  const Point_2 mid(ms / static_cast<int>(k), mt / static_cast<int>(k));

  const std::vector<Linear> E = detail::local_candidates(MB, NB, X, Y, corners, mid, delta, U);
  const std::size_t m = E.size();
  stats.candidates += m;
  stats.largest = std::max(stats.largest, m);

  // First try to show F <= delta via level lines. If F > delta somewhere,
  // delta grows; try again with the new delta (and fewer candidates) a few
  // times, then search for the maximum below.
  if (level_lines > 0 && level_rounds > 0) {
    const LevelOutcome lo = level_certify(corners, E, MB, NB, M, N, delta, level_lines, stats,
                                          evaluated);
    if (lo == LevelOutcome::certified) return LocalOutcome::done;
    if (lo == LevelOutcome::raised)
      return local_polygon(corners, MB, NB, M, N, delta, max_candidates, false, stats, evaluated,
                           level_lines, level_rounds - 1);
  }
  if (m > max_candidates) return LocalOutcome::over_limit;
  std::vector<LinearI> EI;
  for (const Linear& e : E) EI.push_back(LinearI(e));

  // Points where F might exceed delta, with the candidate value there.
  std::vector<std::pair<FT, Point_2> > points;
  auto in_window = [&](const FT& v) { return v > delta && (!bounded || v <= U_exact); };
  const Interval delta_i = IntervalArithmetic::from(delta);
  auto may_window = [&](double sup, double inf) { return sup > delta_i.inf() && inf <= U; };

  // 1. Lines e_i = e_j that cross the face where the common value can lie
  //    in the window. Their boundary points are candidates.
  struct Line {
    std::size_t i, j;
  };
  std::vector<Line> lines;
  std::vector<std::vector<std::size_t> > through(m);  // lines through each function
  for (std::size_t i = 0; i < m; i++)
    for (std::size_t j = i + 1; j < m; j++) {
      const LinearI LI = EI[i] - EI[j];
      std::vector<Interval> v(k);
      bool maybe_pos = false, maybe_neg = false;
      for (std::size_t q = 0; q < k; q++) {
        v[q] = LI.at(si[q], ti[q]);
        maybe_pos = maybe_pos || v[q].sup() > 0;
        maybe_neg = maybe_neg || v[q].inf() < 0;
      }
      if (!maybe_pos || !maybe_neg) continue;
      // Range of e_i over the possible crossing points.
      double hi = -std::numeric_limits<double>::infinity(), lo = std::numeric_limits<double>::infinity();
      bool unsure = false;
      for (std::size_t q = 0; q < k; q++) {
        const std::size_t r = (q + 1) % k;
        if ((v[q].inf() > 0 && v[r].inf() > 0) || (v[q].sup() < 0 && v[r].sup() < 0)) continue;
        const Interval den = v[q] - v[r];
        if (den.inf() <= 0 && den.sup() >= 0) {
          unsure = true;
          continue;
        }
        const Interval w = v[q] / den;
        const Interval val = EI[i].at(si[q] + w * (si[r] - si[q]), ti[q] + w * (ti[r] - ti[q]));
        hi = std::max(hi, val.sup());
        lo = std::min(lo, val.inf());
      }
      if (!unsure && !may_window(hi, lo)) continue;

      // Exact: the segment of the line inside the face.
      const Linear L = E[i] - E[j];
      std::vector<FT> ve(k);
      bool pos = false, neg = false;
      for (std::size_t q = 0; q < k; q++) {
        ve[q] = L.at(corners[q]);
        pos = pos || ve[q] > 0;
        neg = neg || ve[q] < 0;
      }
      if (!pos || !neg) continue;
      std::vector<Point_2> ends;
      for (std::size_t q = 0; q < k; q++) {
        const std::size_t r = (q + 1) % k;
        if (ve[q] == 0) {
          ends.push_back(corners[q]);
        } else if ((ve[q] > 0 && ve[r] < 0) || (ve[q] < 0 && ve[r] > 0)) {
          const FT w = ve[q] / (ve[q] - ve[r]);
          ends.push_back(Point_2(corners[q].x() + w * (corners[r].x() - corners[q].x()),
                                 corners[q].y() + w * (corners[r].y() - corners[q].y())));
        }
      }
      if (ends.size() != 2) continue;  // a convex face meets a crossing line in a segment
      const FT a = E[i].at(ends[0]), b = E[i].at(ends[1]);
      if (!(std::max(a, b) > delta && (!bounded || std::min(a, b) <= U_exact))) continue;
      through[i].push_back(lines.size());
      through[j].push_back(lines.size());
      lines.push_back(Line{i, j});
      if (in_window(a)) points.push_back({a, ends[0]});
      if (in_window(b)) points.push_back({b, ends[1]});
    }
  stats.lines += lines.size();

  // 2. Points inside the face where three candidates e_i = e_j = e_k agree:
  //    intersections of two lines through the same function.
  for (std::size_t i = 0; i < m; i++)
    for (std::size_t x = 0; x < through[i].size(); x++)
      for (std::size_t y = x + 1; y < through[i].size(); y++) {
        const Line& l1 = lines[through[i][x]];
        const Line& l2 = lines[through[i][y]];
        const std::size_t j = l1.i == i ? l1.j : l1.i, h = l2.i == i ? l2.j : l2.i;
        if (j == h) continue;
        const LinearI A1 = EI[i] - EI[j], A2 = EI[i] - EI[h];
        const Interval det = A1.a * A2.b - A1.b * A2.a;
        if (det.inf() > 0 || det.sup() < 0) {
          const Interval ps = (A1.b * A2.c - A2.b * A1.c) / det;
          const Interval pt = (A2.a * A1.c - A1.a * A2.c) / det;
          bool maybe_inside = true;
          for (std::size_t q = 0; q < k && maybe_inside; q++) {
            const std::size_t r = (q + 1) % k;
            const Interval cross = (si[r] - si[q]) * (pt - ti[q]) - (ti[r] - ti[q]) * (ps - si[q]);
            maybe_inside = cross.sup() >= 0;
          }
          const Interval val = EI[i].at(ps, pt);
          if (!maybe_inside || !may_window(val.sup(), val.inf())) continue;
        }
        const Linear B1 = E[i] - E[j], B2 = E[i] - E[h];
        const FT d = B1.a * B2.b - B1.b * B2.a;
        if (d == 0) continue;  // parallel or the same line: covered by its ends
        const Point_2 p((B1.b * B2.c - B2.b * B1.c) / d, (B2.a * B1.c - B1.a * B2.c) / d);
        bool inside = true;
        for (std::size_t q = 0; q < k && inside; q++) {
          const std::size_t r = (q + 1) % k;
          inside = CGAL::orientation(corners[q], corners[r], p) != CGAL::RIGHT_TURN;
        }
        if (!inside) continue;
        const FT val = E[i].at(p);
        if (in_window(val)) points.push_back({val, p});
      }
  stats.points += points.size();

  // 3. Evaluate, best candidate value first. F at a point can only exceed
  //    delta if its candidate value does. At a single point the bar-box
  //    boxes are the bars themselves, so its test is a cheap exact check
  //    of F <= delta; only points it cannot certify are evaluated.
  std::sort(points.begin(), points.end(),
            [](const std::pair<FT, Point_2>& x, const std::pair<FT, Point_2>& y) {
              return x.first > y.first;
            });
  for (const auto& vp : points) {
    if (!(vp.first > delta)) break;
    evaluate(vp.second);
  }
  return LocalOutcome::done;
}

namespace detail {

// The convex polygon cut by the half-plane a s + b t + c >= 0 (exact).
inline std::vector<Point_2> clip(const std::vector<Point_2>& poly, const FT& a, const FT& b,
                                 const FT& c) {
  std::vector<Point_2> out;
  const std::size_t k = poly.size();
  for (std::size_t q = 0; q < k; q++) {
    const Point_2& p0 = poly[q];
    const Point_2& p1 = poly[(q + 1) % k];
    const FT v0 = a * p0.x() + b * p0.y() + c, v1 = a * p1.x() + b * p1.y() + c;
    if (v0 >= 0) out.push_back(p0);
    if ((v0 > 0 && v1 < 0) || (v0 < 0 && v1 > 0)) {
      const FT w = v0 / (v0 - v1);
      out.push_back(Point_2(p0.x() + w * (p1.x() - p0.x()), p0.y() + w * (p1.y() - p0.y())));
    }
  }
  // Drop repeated points (a cut through a corner).
  std::vector<Point_2> clean;
  for (const Point_2& p : out)
    if (clean.empty() || p != clean.back()) clean.push_back(p);
  while (clean.size() > 1 && clean.front() == clean.back()) clean.pop_back();
  return clean;
}

}  // namespace detail

// The local phase on a polygon, splitting it in two while it has too many
// candidates, at most max_depth times. A piece of a face of A0 lies in that
// face, so the argument of the local phase applies to it unchanged; only its
// new corners need checking. The cut halves the larger of the two extents
// X (s_max - s_min) and t_max - t_min, which bound how far pushes move
// (X = largest |v_x| of the strip), so the push ranges of the pieces shrink.
inline void local_region(const std::vector<Point_2>& corners, const BarBoxModule& MB,
                         const BarBoxModule& NB, const Module& M, const Module& N,
                         const FT& max_x, FT& delta, std::size_t max_candidates,
                         std::size_t depth, std::size_t max_depth, bool new_corners,
                         LocalStats& stats, std::set<Point_2, PointLess>& evaluated,
                         std::size_t level_lines = 0) {
  const LocalOutcome outcome = local_polygon(corners, MB, NB, M, N, delta, max_candidates,
                                             new_corners, stats, evaluated, level_lines);
  if (outcome != LocalOutcome::over_limit) return;
  if (depth >= max_depth) {
    ++stats.over_limit;
    return;
  }
  ++stats.splits;
  FT smin = corners[0].x(), smax = smin, tmin = corners[0].y(), tmax = tmin;
  for (const Point_2& c : corners) {
    smin = std::min(smin, c.x());
    smax = std::max(smax, c.x());
    tmin = std::min(tmin, c.y());
    tmax = std::max(tmax, c.y());
  }
  // Cut s = s_mid or t = t_mid, as halves of the polygon.
  const bool cut_s = max_x * (smax - smin) >= tmax - tmin;
  const FT mid = cut_s ? (smin + smax) / 2 : (tmin + tmax) / 2;
  const FT a = cut_s ? FT(1) : FT(0), b = cut_s ? FT(0) : FT(1);
  for (int side = 0; side < 2; side++) {
    const FT sign = side ? FT(-1) : FT(1);
    const std::vector<Point_2> piece = detail::clip(corners, sign * a, sign * b, -sign * mid);
    if (piece.size() < 3) continue;
    local_region(piece, MB, NB, M, N, max_x, delta, max_candidates, depth + 1, max_depth, true,
                 stats, evaluated, level_lines);
  }
}

// Examines the red faces of a strip's colored arrangement (after
// insert_strip) and raises delta to the maximum of F over them. M and N are
// the strip's modules (swapped for the steep strip). Faces with more than
// max_candidates candidates are split up to max_depth times (0: never).
inline void local_phase(const Colored_arrangement& arr, const Module& M, const Module& N,
                        FT& delta, std::size_t max_candidates, LocalStats& stats,
                        std::size_t max_depth = 0, std::size_t level_lines = 0) {
  const BarBoxModule MB = bar_box_module(M), NB = bar_box_module(N);
  const FT X = max_abs_x(M, N);
  std::set<Point_2, PointLess> evaluated;
  for (auto f = arr.faces_begin(); f != arr.faces_end(); ++f) {
    if (f->is_unbounded() || f->data() != Color::red) continue;
    const std::vector<Point_2> corners = face_corners(f);
    if (bar_box_green(MB, NB, corners, delta)) continue;  // delta may have grown
    ++stats.red_faces;
    local_region(corners, MB, NB, M, N, X, delta, max_candidates, 0, max_depth, false, stats,
                 evaluated, level_lines);
  }
}

// --- The whole algorithm ------------------------------------------------------

struct IncrementalOptions {
  std::size_t samples = 50;            // random A0 vertices evaluated first (half per strip, at least one)
  unsigned seed = 1;                   // of the sampling and of the insertion order
  bool skip_covered = true;            // do not evaluate vertices covered by a sample
  FaceBound bound = FaceBound::both;   // face bound of the insertion
  bool local_phase = false;            // go beyond the A0 vertices: the matching distance
  std::size_t max_candidates = 500;    // local phase: skip faces with more candidates
  std::size_t split_depth = 0;         // local phase: split such faces up to this depth first
  std::size_t level_lines = 0;         // local phase: certify via at most this many level lines (0: off)
};

struct IncrementalResult {
  // Maximum of F over the A0 vertices; std::nullopt means infinity.
  std::optional<FT> a0_max;
  // With local_phase: the matching distance, or only a lower bound if
  // local.over_limit > 0. Equal to a0_max otherwise.
  std::optional<FT> distance;
  StripStats flat, steep;
  LocalStats local;
};

// The maximum of F over the A0 vertices of both strips, computed with the
// red/green insertion, and with options.local_phase the matching distance.
// Infinite (std::nullopt) if the modules have different numbers of infinite
// bars. First evaluates F at random A0 vertices for a good starting Delta,
// then inserts the segments of the flat and then the steep strip; Delta is
// shared by both. The A0 result equals max_F_on_vertices over both strips.
inline IncrementalResult incremental_max_F(const Module& M, const Module& N,
                                           const IncrementalOptions& options = IncrementalOptions()) {
  IncrementalResult result;
  const Module M_swap = swapped_module(M), N_swap = swapped_module(N);
  std::mt19937 rng(options.seed);
  FT delta(0);
  // The first sample also decides whether F is infinite (everywhere).
  std::vector<Sample> flat_samples, steep_samples;
  if (!sample_delta(M, N, std::max<std::size_t>(1, (options.samples + 1) / 2), rng, delta,
                    flat_samples, result.flat))
    return result;
  sample_delta(M_swap, N_swap, options.samples / 2, rng, delta, steep_samples, result.steep);

  Colored_arrangement flat, steep;
  insert_strip(flat, M, N, rng, delta, result.flat,
               options.skip_covered ? &flat_samples : nullptr, options.bound);
  insert_strip(steep, M_swap, N_swap, rng, delta, result.steep,
               options.skip_covered ? &steep_samples : nullptr, options.bound);
  result.a0_max = delta;

  if (options.local_phase) {
    local_phase(flat, M, N, delta, options.max_candidates, result.local, options.split_depth,
                options.level_lines);
    local_phase(steep, M_swap, N_swap, delta, options.max_candidates, result.local,
                options.split_depth, options.level_lines);
  }
  result.distance = delta;
  return result;
}

}  // namespace matching
