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
#include <CGAL/Arr_observer.h>

#include <random>

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
  Coloring_observer(Colored_arrangement& arr, const Module& M, const Module& N, FT& delta,
                    const std::vector<Sample>* samples = nullptr)
      : CGAL::Arr_observer<Colored_arrangement>(arr), M_(M), N_(N), delta_(delta),
        max_x_(max_abs_x(M, N)), samples_(samples) {}

  // X of the face bound for this strip (see max_abs_x).
  const FT& max_x() const { return max_x_; }

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
};

inline FT face_upper_bound(Colored_arrangement::Face_handle f, const FT& max_x) {
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

// Inserts one segment, then restores the invariant: no yellow faces, and
// no edge with green on both sides.
inline void insert_segment(Colored_arrangement& arr, Coloring_observer& obs,
                           const Segment_2& seg, const FT& delta) {
  CGAL::insert(arr, X_monotone_curve_2(seg));   // the callbacks fill obs.yellow and obs.green_edges

  // Step 1: classify every yellow face as red or green (face_upper_bound vs delta).
  std::vector<Colored_arrangement::Face_handle> new_green_faces;
  for (Colored_arrangement::Face_handle f : obs.yellow)
    if (face_upper_bound(f, obs.max_x()) <= delta) {
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
// Coloring_observer).
inline void insert_strip(const Module& M, const Module& N, std::mt19937& rng, FT& delta,
                         StripStats& stats, const std::vector<Sample>* samples = nullptr) {
  Colored_arrangement arr;
  build_initial_quadrilateral(arr, M, N, delta);
  stats.evaluations = arr.number_of_vertices();  // the corners

  // Declared after arr, so it is destroyed (and detached) before arr.
  Coloring_observer obs(arr, M, N, delta, samples);

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

// The maximum of F over the A0 vertices of both strips, computed with the
// red/green insertion. Returns std::nullopt (infinity) if the modules have
// different numbers of infinite bars. First evaluates F at `samples` random
// A0 vertices (half per strip, at least one) for a good starting Delta,
// then inserts the segments of the flat and then the steep strip. Delta is
// shared by both. With skip_covered, vertices covered by a sample are not
// evaluated. The result equals max_F_on_vertices over both strips.
inline std::optional<FT> incremental_max_F(const Module& M, const Module& N,
                                           std::size_t samples, unsigned seed,
                                           StripStats& flat, StripStats& steep,
                                           bool skip_covered = true) {
  const Module M_swap = swapped_module(M), N_swap = swapped_module(N);
  std::mt19937 rng(seed);
  FT delta(0);
  // The first sample also decides whether F is infinite (everywhere).
  std::vector<Sample> flat_samples, steep_samples;
  if (!sample_delta(M, N, std::max<std::size_t>(1, (samples + 1) / 2), rng, delta,
                    flat_samples, flat))
    return std::nullopt;
  sample_delta(M_swap, N_swap, samples / 2, rng, delta, steep_samples, steep);
  insert_strip(M, N, rng, delta, flat, skip_covered ? &flat_samples : nullptr);
  insert_strip(M_swap, N_swap, rng, delta, steep, skip_covered ? &steep_samples : nullptr);
  return delta;
}

}  // namespace matching
