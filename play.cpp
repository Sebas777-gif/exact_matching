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
              << " [--incremental [--samples k] [--seed n] [--no-skip]]"
              << " (each file .scc, .scc_k, .minpres_k, .sccsum, or .bif)"
              << std::endl;
    std::exit(1);
  }

  std::vector<std::string> input_files;
  bool preprocess = true;
  bool evaluate = true;
  bool incremental = false;
  bool skip_covered = true;
  long samples = 50;
  long seed = 1;
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
    } else if (arg == "--samples" || arg == "--seed") {
      if (i + 1 >= argc) {
        std::cerr << arg << " needs an integer" << std::endl;
        std::exit(1);
      }
      const long value = std::stol(argv[++i]);
      if (value < 0) {
        std::cerr << arg << " must be non-negative" << std::endl;
        std::exit(1);
      }
      (arg == "--samples" ? samples : seed) = value;
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
    matching::StripStats flat_stats, steep_stats;
    const std::optional<matching::FT> result = matching::incremental_max_F(
        modules[0], modules[1], samples, seed, flat_stats, steep_stats, skip_covered);
    auto strip = [](const char* name, const matching::StripStats& st) {
      std::cout << name << st.samples << " samples, " << st.evaluations << " evaluations, "
                << st.skipped << " skipped, " << st.segments << " segments; final " << st.vertices << " vertices, "
                << st.edges << " edges, " << st.faces << " faces (" << st.red_faces
                << " red)" << std::endl;
    };
    std::cout << "Incremental: " << std::max(1L, samples) << " samples, seed " << seed
              << (skip_covered ? ", skipping covered vertices" : "") << std::endl;
    if (result) {
      strip("Flat strip:  ", flat_stats);
      strip("Steep strip: ", steep_stats);
    }
    const std::size_t evaluations = flat_stats.samples + flat_stats.evaluations
        + steep_stats.samples + steep_stats.evaluations;
    std::cout << "F evaluations: " << evaluations << " (A0 has " << a0_vertices
              << " vertices)" << std::endl;
    std::cout << "Max F on A0 vertices: " << format(result) << std::endl;
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
