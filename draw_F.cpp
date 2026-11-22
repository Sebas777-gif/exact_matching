// Draws the A0 arrangement of a strip with CGAL's basic viewer (Qt6), each
// bounded face colored by F at its centroid: a light-to-dark blue ramp from
// the smallest face value up to Delta, and orange where F equals Delta, the
// maximum of F over the A0 vertices of both strips (exact comparison). So the
// orange faces are the plateaus on which F = Delta, where bounds of the form
// "F at a corner plus movement" cannot turn a face green.
//
// The viewer has no legend; draw_F prints the color scale to stdout. The
// walls s = 0 and s = 1 are labeled at both ends, with t ticks along them
// (for the steep strip, s and t are those of the swapped modules: slices
// x = s y + t). In the window, the mouse pans and zooms; T toggles the
// labels, H shows the viewer's key bindings.
// F is evaluated once per vertex (for Delta) and once per bounded face.
//
// Usage: draw_F file1 file2 [--dim d] [--no-aida] [--strip flat|steep|both] [--no-scale]
//               [--no-window]
//   --no-window  only print the statistics (plateau faces and their share of the area)
//   --strip   which strip to show (default both: flat first, steep after the
//             flat window is closed)
//   --no-scale  draw (s, t) as is. By default s is stretched by X, the largest
//             |v_x| of the strip, so that a unit step in s moves a push about
//             as far as a unit step in t (otherwise the strip of width 1 is
//             a thin sliver next to t ranges in the tens or hundreds).

#include <CGAL/draw_arrangement_2.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <aida_preprocess.hpp>
#include <matching_distance.h>

namespace {

using matching::Arrangement_2;
using matching::FT;
using matching::Module;
using matching::Point_2;

// Sequential blue ramp, steps 100..700 (light to dark).
const char* const kBlue[] = {"#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec",
                             "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab",
                             "#184f95", "#104281", "#0d366b"};
const std::size_t kSteps = sizeof(kBlue) / sizeof(kBlue[0]);
const char* const kPlateau = "#eb6834";  // orange: F = Delta

CGAL::IO::Color hex(const char* h) {
  auto byte = [&](int i) { return static_cast<unsigned char>(std::stoi(std::string(h + i, 2), nullptr, 16)); };
  return CGAL::IO::Color(byte(1), byte(3), byte(5));
}

// The A0 arrangement of a strip with s stretched by `scale` (an affine map,
// so the faces are the same as in build_arrangement_A).
void build_scaled(Arrangement_2& arr, const Module& M, const Module& N, const FT& scale) {
  std::set<matching::Grade2> grades;
  matching::collect_grades(M, N, grades);
  std::vector<matching::X_monotone_curve_2> curves;
  FT t0_min = grades.begin()->y, t0_max = t0_min;
  FT t1_min = grades.begin()->y - grades.begin()->x, t1_max = t1_min;
  for (const matching::Grade2& g : grades) {
    const FT t0 = g.y, t1 = g.y - g.x;
    curves.emplace_back(matching::Segment_2(Point_2(FT(0), t0), Point_2(scale, t1)));
    t0_min = std::min(t0_min, t0);
    t0_max = std::max(t0_max, t0);
    t1_min = std::min(t1_min, t1);
    t1_max = std::max(t1_max, t1);
  }
  if (t0_min < t0_max) curves.emplace_back(matching::Segment_2(Point_2(FT(0), t0_min), Point_2(FT(0), t0_max)));
  if (t1_min < t1_max) curves.emplace_back(matching::Segment_2(Point_2(scale, t1_min), Point_2(scale, t1_max)));
  CGAL::insert(arr, curves.begin(), curves.end());
}

FT max_abs_x(const Module& M, const Module& N) {
  std::set<matching::Grade2> grades;
  matching::collect_grades(M, N, grades);
  FT x(0);
  for (const matching::Grade2& g : grades) x = std::max(x, CGAL::abs(g.x));
  return x;
}

// Tick values for [lo, hi]: multiples of a step 1, 2 or 5 times a power of
// ten, about `count` of them.
std::vector<double> nice_ticks(double lo, double hi, int count) {
  std::vector<double> ticks;
  if (!(hi > lo)) {
    ticks.push_back(lo);
    return ticks;
  }
  const double raw = (hi - lo) / count, mag = std::pow(10.0, std::floor(std::log10(raw)));
  const double step = raw / mag < 1.5 ? mag : raw / mag < 3.5 ? 2 * mag : raw / mag < 7.5 ? 5 * mag : 10 * mag;
  for (double t = std::ceil(lo / step) * step; t <= hi + 1e-9 * step; t += step)
    ticks.push_back(std::abs(t) < 1e-9 * step ? 0.0 : t);
  return ticks;
}

std::string number(double x) {
  std::ostringstream out;
  out << x;
  return out.str();
}

// Orientation labels: "s = 0" and "s = 1" at both ends of the two walls, and
// t ticks with values along them (left of s = 0, right of s = 1). The viewer
// draws text in screen space, so it stays readable when zooming; T toggles it.
template <class Arr>
void add_labels(CGAL::Graphics_scene& scene, const Arr& arr, double scale) {
  double t0_lo = 1e300, t0_hi = -1e300, t1_lo = 1e300, t1_hi = -1e300;
  for (auto v = arr.vertices_begin(); v != arr.vertices_end(); ++v) {
    const double x = CGAL::to_double(v->point().x()), t = CGAL::to_double(v->point().y());
    if (x == 0) {
      t0_lo = std::min(t0_lo, t);
      t0_hi = std::max(t0_hi, t);
    }
    if (x == scale) {
      t1_lo = std::min(t1_lo, t);
      t1_hi = std::max(t1_hi, t);
    }
  }
  if (t0_lo > t0_hi || t1_lo > t1_hi) return;
  const double span = std::max({t0_hi, t1_hi}) - std::min({t0_lo, t1_lo});
  const double gap = 0.04 * std::max(span, scale), tick = 0.015 * std::max(span, scale);
  typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_2 P;
  const CGAL::IO::Color ink(60, 60, 60);
  scene.add_text(P(0, t0_lo - gap), "s = 0");
  scene.add_text(P(0, t0_hi + gap), "s = 0");
  scene.add_text(P(scale, t1_lo - gap), "s = 1");
  scene.add_text(P(scale, t1_hi + gap), "s = 1");
  const double lo = std::min(t0_lo, t1_lo), hi = std::max(t0_hi, t1_hi);
  for (double t : nice_ticks(lo, hi, 8)) {
    if (t >= t0_lo && t <= t0_hi) {
      scene.add_segment(P(-tick, t), P(0, t), ink);
      // Text starts at its anchor, so a label left of the wall needs room
      // for its own width (about two gaps at the initial zoom).
      scene.add_text(P(-tick - 2.4 * gap, t), "t = " + number(t));
    }
    if (t >= t1_lo && t <= t1_hi) {
      scene.add_segment(P(scale, t), P(scale + tick, t), ink);
      scene.add_text(P(scale + tick + gap, t), "t = " + number(t));
    }
  }
}

struct Strip {
  const char* name;
  Module M, N;
  FT scale;
  Arrangement_2 arr;
  std::map<Arrangement_2::Face_const_handle, FT> value;  // F at the centroid
};

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> files;
  bool preprocess = true, scaled = true, window = true;
  int dim = 0;
  std::string which = "both";
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--no-aida") preprocess = false;
    else if (arg == "--no-scale") scaled = false;
    else if (arg == "--no-window") window = false;
    else if (arg == "--dim" && i + 1 < argc) dim = std::stoi(argv[++i]);
    else if (arg == "--strip" && i + 1 < argc) which = argv[++i];
    else if (arg.rfind("--", 0) == 0) {
      std::cerr << "Unknown flag " << arg << std::endl;
      return 1;
    } else files.push_back(arg);
  }
  if (files.size() != 2 || (which != "flat" && which != "steep" && which != "both")) {
    std::cerr << "Usage: " << argv[0]
              << " file1 file2 [--dim d] [--no-aida] [--strip flat|steep|both] [--no-scale] [--no-window]"
              << std::endl;
    return 1;
  }
  Module M, N;
  try {
    M = scc_hera::load_module(files[0], dim, preprocess);
    N = scc_hera::load_module(files[1], dim, preprocess);
  } catch (const std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    return 1;
  }

  std::vector<Strip> strips(2);
  strips[0].name = "flat";
  strips[0].M = M;
  strips[0].N = N;
  strips[1].name = "steep";
  strips[1].M = matching::swapped_module(M);
  strips[1].N = matching::swapped_module(N);

  // Delta over the vertices of both strips, then F at every face centroid.
  FT delta(0);
  for (Strip& st : strips) {
    st.scale = scaled ? std::max(FT(1), max_abs_x(st.M, st.N)) : FT(1);
    build_scaled(st.arr, st.M, st.N, st.scale);
    for (auto v = st.arr.vertices_begin(); v != st.arr.vertices_end(); ++v) {
      const std::optional<FT> f = matching::F(st.M, st.N, v->point().x() / st.scale, v->point().y());
      if (!f) {
        std::cout << "F is infinite (different numbers of infinite bars); nothing to draw." << std::endl;
        return 0;
      }
      delta = std::max(delta, *f);
    }
  }
  FT lo = delta;
  for (Strip& st : strips) {
    for (auto f = st.arr.faces_begin(); f != st.arr.faces_end(); ++f) {
      if (f->is_unbounded()) continue;
      FT ms(0), mt(0);
      int n = 0;
      auto h = f->outer_ccb(), start = h;
      do {
        ms += h->source()->point().x();
        mt += h->source()->point().y();
        ++n;
      } while (++h != start);
      // Faces of A0 are convex, so the average of the corners is inside.
      const FT value = *matching::F(st.M, st.N, ms / n / st.scale, mt / n);
      st.value[f] = value;
      lo = std::min(lo, value);
    }
  }

  // Color scale, printed since the viewer has no legend.
  const double dlo = CGAL::to_double(lo), dhi = CGAL::to_double(delta);
  std::cout << "Delta (max F over the A0 vertices of both strips): " << CGAL::exact(delta) << " ("
            << dhi << ")" << std::endl;
  std::cout << "Face color = F at the face centroid. Orange " << kPlateau << ": F = Delta (plateau). "
            << "Blue ramp, light to dark: " << dlo << " to " << dhi << "." << std::endl;
  for (const Strip& st : strips) {
    // Plateau share by area, in the (s, t) plane (the drawing stretches s).
    std::size_t plateau = 0;
    double area = 0, plateau_area = 0;
    for (const auto& kv : st.value) {
      double a = 0;
      auto h = kv.first->outer_ccb(), start = h;
      do {
        const Point_2& p = h->source()->point();
        const Point_2& q = h->target()->point();
        a += CGAL::to_double(p.x() * q.y() - q.x() * p.y());
      } while (++h != start);
      a = std::abs(a) / 2 / CGAL::to_double(st.scale);
      area += a;
      if (kv.second == delta) {
        ++plateau;
        plateau_area += a;
      }
    }
    std::cout << "  " << std::setw(5) << st.name << " strip: " << st.value.size() << " faces, "
              << plateau << " on the plateau F = Delta (" << std::fixed << std::setprecision(1)
              << (area > 0 ? 100 * plateau_area / area : 0.0) << "% of the area)"
              << std::defaultfloat << std::setprecision(6)
              << (scaled ? ", s drawn stretched by X = " + number(CGAL::to_double(st.scale)) : std::string())
              << std::endl;
  }
  if (!window) return 0;

  for (Strip& st : strips) {
    if (which != "both" && which != st.name) continue;
    CGAL::Graphics_scene_options<Arrangement_2, Arrangement_2::Vertex_const_handle,
                                 Arrangement_2::Halfedge_const_handle,
                                 Arrangement_2::Face_const_handle> gso;
    gso.colored_face = [](const Arrangement_2&, Arrangement_2::Face_const_handle) { return true; };
    gso.face_color = [&](const Arrangement_2&, Arrangement_2::Face_const_handle f) {
      auto it = st.value.find(f);
      if (it == st.value.end()) return hex(kBlue[0]);
      if (it->second == delta) return hex(kPlateau);
      const double x = dhi > dlo ? (CGAL::to_double(it->second) - dlo) / (dhi - dlo) : 0.0;
      const std::size_t k = std::min(kSteps - 1, static_cast<std::size_t>(x * (kSteps - 1) + 0.5));
      return hex(kBlue[k]);
    };
    const std::string title = std::string("F on the ") + st.name + " strip"
        + (std::string(st.name) == "flat" ? " (slices y = s x + t)" : " (slices x = s y + t)")
        + ", orange: F = Delta, blue: lighter = smaller F";
    // As CGAL::draw, but with a larger window than the viewer's default.
    CGAL::Graphics_scene scene;
    CGAL::add_to_graphics_scene(st.arr, scene, gso);
    add_labels(scene, st.arr, CGAL::to_double(st.scale));
    CGAL::Qt::init_ogl_context(4, 3);
    int qt_argc = 1;
    const char* qt_argv[2] = {title.c_str(), nullptr};
    QApplication app(qt_argc, const_cast<char**>(qt_argv));
    CGAL::Qt::Basic_viewer viewer(app.activeWindow(), scene, title.c_str());
    viewer.resize(1100, 850);
    viewer.show();
    app.exec();
  }
  return 0;
}
