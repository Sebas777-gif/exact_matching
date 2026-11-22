#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <colored_arrangement.h>
#include <matching_distance.h>
#include <aida_preprocess.hpp>

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " file1 file2 [--dim d] [--no-aida] [--a0-only]"
              << " [--incremental [--samples k] [--seed n] [--no-skip]"
              << " [--bound lipschitz|variation|barbox|both]]"
              << " [--distance [--max-candidates k] [--split-depth d] [--level-lines l]]"
              << " (each file .scc, .scc_k, .minpres_k, .sccsum, or .bif)"
              << std::endl;
    std::exit(1);
  }

  std::vector<std::string> input_files;
  bool preprocess = true;
  bool evaluate = true;
  bool incremental = false;
  bool skip_covered = true;
  matching::FaceBound bound = matching::FaceBound::both;
  bool distance = false;
  long samples = 50;
  long seed = 1;
  long max_candidates = 500;
  long split_depth = 0;
  long level_lines = 0;
  int homology_dim = 0;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--no-aida") {
      preprocess = false;
    } else if (arg == "--a0-only") {
      evaluate = false;
    } else if (arg == "--incremental") {
      incremental = true;
    } else if (arg == "--no-skip") {
      skip_covered = false;
    } else if (arg == "--bound") {
      const std::string kind = i + 1 < argc ? argv[++i] : "";
      if (kind == "lipschitz") {
        bound = matching::FaceBound::lipschitz;
      } else if (kind == "variation") {
        bound = matching::FaceBound::variation;
      } else if (kind == "barbox") {
        bound = matching::FaceBound::barbox;
      } else if (kind == "both") {
        bound = matching::FaceBound::both;
      } else {
        std::cerr << "--bound needs lipschitz, variation, barbox, or both" << std::endl;
        std::exit(1);
      }
    } else if (arg == "--distance") {
      distance = true;
      incremental = true;  // the local phase works on the colored arrangement
    } else if (arg == "--samples" || arg == "--seed" || arg == "--max-candidates"
               || arg == "--split-depth" || arg == "--level-lines") {
      if (i + 1 >= argc) {
        std::cerr << arg << " needs an integer" << std::endl;
        std::exit(1);
      }
      const long value = std::stol(argv[++i]);
      if (value < 0) {
        std::cerr << arg << " must be non-negative" << std::endl;
        std::exit(1);
      }
      (arg == "--samples" ? samples : arg == "--seed" ? seed
       : arg == "--max-candidates" ? max_candidates
       : arg == "--split-depth" ? split_depth : level_lines) = value;
    } else if (arg == "--dim") {
      if (i + 1 >= argc) {
        std::cerr << "--dim needs an integer" << std::endl;
        std::exit(1);
      }
      homology_dim = std::stoi(argv[++i]);
      if (homology_dim < 0) {
        std::cerr << "--dim must be non-negative" << std::endl;
        std::exit(1);
      }
    } else if (arg.size() >= 2 && arg[0] == '-' && arg[1] == '-') {
      std::cerr << "Unknown flag " << arg << std::endl;
      std::exit(1);
    } else {
      input_files.push_back(arg);
    }
  }

  if (input_files.size() != 2) {
    std::cerr << "Expected two input files (one module each), got "
              << input_files.size() << std::endl;
    std::exit(1);
  }

  std::vector<matching::Module> modules;
  try {
    for (const std::string& input_filename : input_files) {
      modules.push_back(scc_hera::load_module(input_filename, homology_dim, preprocess));
    }
  } catch (const std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    std::exit(1);
  }

  std::cout << "Loaded 2 module(s) with " << modules[0].summands.size()
            << " summand(s) and " << modules[1].summands.size() << " summand(s)"
            << std::endl;

  // Flat strip: slices y = s x + t with 0 <= s <= 1.
  // Steep strip: slices x = s y + t with 0 <= s <= 1, i.e. the flat strip
  // of the modules with x and y swapped.
  matching::Arrangement_2 flat, steep;
  const std::size_t nlines = matching::build_arrangement_A(flat, modules[0], modules[1]);
  matching::build_arrangement_A(steep, matching::swapped_module(modules[0]),
                                matching::swapped_module(modules[1]));

  auto counts = [](const char* name, const matching::Arrangement_2& arr) {
    std::cout << name << arr.number_of_vertices() << " vertices, "
              << arr.number_of_edges() << " edges, "
              << arr.number_of_faces() << " faces" << std::endl;
  };
  std::cout << "A0 lines: " << nlines << std::endl;
  counts("Flat strip:  ", flat);
  counts("Steep strip: ", steep);
  std::cout << "Arrangement has "
            << flat.number_of_vertices() + steep.number_of_vertices() << " vertices, "
            << flat.number_of_edges() + steep.number_of_edges() << " edges, and "
            << flat.number_of_faces() + steep.number_of_faces() << " faces" << std::endl;
  if (!evaluate) return 0;

  // Exact value as a fraction ("3/2"), then a decimal; "inf" for infinity.
  auto format = [](const std::optional<matching::FT>& v) -> std::string {
    if (!v) return "inf";
    std::ostringstream out;
    out << CGAL::exact(*v) << " (" << CGAL::to_double(*v) << ")";
    return out.str();
  };
  const std::size_t a0_vertices = flat.number_of_vertices() + steep.number_of_vertices();

  if (incremental) {
    matching::IncrementalOptions options;
    options.samples = samples;
    options.seed = seed;
    options.skip_covered = skip_covered;
    options.bound = bound;
    options.local_phase = distance;
    options.max_candidates = max_candidates;
    options.split_depth = split_depth;
    options.level_lines = level_lines;
    const matching::IncrementalResult result =
        matching::incremental_max_F(modules[0], modules[1], options);

    auto strip = [](const char* name, const matching::StripStats& st) {
      std::cout << name << st.samples << " samples, " << st.evaluations << " evaluations, "
                << st.skipped << " skipped, " << st.segments << " segments; final " << st.vertices << " vertices, "
                << st.edges << " edges, " << st.faces << " faces (" << st.red_faces
                << " red)" << std::endl;
    };
    const char* bound_name = bound == matching::FaceBound::lipschitz ? "Lipschitz"
                           : bound == matching::FaceBound::variation ? "variation"
                           : bound == matching::FaceBound::barbox    ? "bar-box"
                                                                     : "variation and bar-box";
    std::cout << "Incremental: " << std::max(1L, samples) << " samples, seed " << seed
              << (skip_covered ? ", skipping covered vertices" : "") << ", " << bound_name
              << " bound" << std::endl;
    if (result.a0_max) {
      strip("Flat strip:  ", result.flat);
      strip("Steep strip: ", result.steep);
    }
    const std::size_t evaluations = result.flat.samples + result.flat.evaluations
        + result.steep.samples + result.steep.evaluations;
    std::cout << "F evaluations: " << evaluations << " (A0 has " << a0_vertices
              << " vertices)" << std::endl;
    // Evaluating F at every A0 vertex is the brute force; this is what the
    // pruning saved (negative if the samples cost more than they saved).
    const double pruned = a0_vertices ? 100.0 * (1.0 - double(evaluations) / a0_vertices) : 0.0;
    std::cout << "Pruned: " << std::fixed << std::setprecision(1) << pruned << std::defaultfloat
              << std::setprecision(6) << "% of the A0 vertex evaluations avoided" << std::endl;
    std::cout << "Max F on A0 vertices: " << format(result.a0_max) << std::endl;
    if (distance) {
      const matching::LocalStats& L = result.local;
      std::cout << "Local phase: " << L.red_faces << " red faces, " << L.candidates
                << " candidate functions (at most " << L.largest << " in a face), " << L.lines
                << " relevant lines, " << L.points << " candidate points, " << L.certified
                << " certified without evaluating, " << L.evaluations
                << " evaluations, " << L.splits << " splits, " << L.over_limit
                << " faces over the candidate limit"
                << std::endl;
      if (level_lines)
        std::cout << "Level lines: " << L.level_tries << " polygons tried, " << L.level_certified
                  << " certified, " << L.level_raised << " raised Delta, " << L.level_skipped
                  << " over the line limit; " << L.level_lines << " lines, " << L.level_cells
                  << " cells, " << L.level_checks << " cells checked at a point" << std::endl;
      std::cout << "Matching distance: " << format(result.distance);
      if (L.over_limit)
        std::cout << ", lower bound (" << L.over_limit << " faces over the candidate limit)";
      else
        std::cout << ", exact";
      std::cout << std::endl;
    }
    return 0;
  }

  const std::optional<matching::FT> max_flat =
      matching::max_F_on_vertices(flat, modules[0], modules[1]);
  const std::optional<matching::FT> max_steep = matching::max_F_on_vertices(
      steep, matching::swapped_module(modules[0]), matching::swapped_module(modules[1]));
  std::optional<matching::FT> max_both;
  if (max_flat && max_steep) max_both = std::max(*max_flat, *max_steep);
  std::cout << "Max F on flat strip vertices:  " << format(max_flat) << std::endl;
  std::cout << "Max F on steep strip vertices: " << format(max_steep) << std::endl;
  // The brute force evaluates every vertex once per strip and stops at the
  // first infinite value.
  const std::size_t evaluations = (max_flat ? flat.number_of_vertices() : 1)
      + (max_steep ? steep.number_of_vertices() : 1);
  std::cout << "F evaluations: " << evaluations << " (A0 has "
            << a0_vertices << " vertices)" << std::endl;
  std::cout << "Max F on A0 vertices: " << format(max_both) << std::endl;
  return 0;
}
