#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <matching_distance.h>
#include <aida_preprocess.hpp>

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " file1 file2 [--dim d] [--no-aida]"
              << " (each file .scc, .scc_k, .minpres_k, .sccsum, or .bif)"
              << std::endl;
    std::exit(1);
  }

  std::vector<std::string> input_files;
  bool preprocess = true;
  int homology_dim = 0;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--no-aida") {
      preprocess = false;
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
  return 0;
}
