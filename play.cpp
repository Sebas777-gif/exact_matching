#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <matching_distance.h>
#include <aida_preprocess.hpp>

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " file.scc|file.scc_k|file.sccsum|file.bif [file2] [--dim d] [--no-aida]"
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

  if (input_files.empty()) {
    std::cerr << "No input file given" << std::endl;
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

  // One file with exactly two presentations: the two scc2020 blocks are the two modules.
  if (modules.size() == 1 && modules[0].summands.size() == 2) {
    matching::Module A = matching::as_module(modules[0].summands[0]);
    matching::Module B = matching::as_module(modules[0].summands[1]);
    modules.clear();
    modules.push_back(std::move(A));
    modules.push_back(std::move(B));
  }

  std::cout << "Loaded " << modules.size() << " module(s)";
  if (!modules.empty()) {
    std::cout << " with";
    for (std::size_t i = 0; i < modules.size(); ++i) {
      if (i) std::cout << " and";
      std::cout << " " << modules[i].summands.size() << " summand(s)";
    }
  }
  std::cout << std::endl;

  if (modules.size() != 2) {
    std::cerr << "A0 needs two modules (two files, or one file with two scc2020 blocks)"
              << std::endl;
    std::exit(1);
  }

  matching::Arrangement_2 arr;
  const std::size_t nlines = matching::build_arrangement_A(arr, modules[0], modules[1]);
  std::cout << "A0 lines: " << nlines << std::endl;
  std::cout << "Arrangement has " << arr.number_of_vertices() << " vertices, "
            << arr.number_of_edges() << " edges, and "
            << arr.number_of_faces() << " faces" << std::endl;
  return 0;
}
