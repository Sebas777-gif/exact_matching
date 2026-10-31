#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <matching_distance.h>
#include <aida_preprocess.hpp>

typedef matching::Kernel Kernel;
typedef matching::Traits_2 Traits_2;
typedef matching::Point_2 Point_2;
typedef matching::Ray_2 Ray_2;
typedef matching::Line_2 Line_2;
typedef matching::X_monotone_curve_2 X_monotone_curve_2;
typedef matching::Arrangement_2 Arrangement_2;
typedef Arrangement_2::Vertex_handle Vertex_handle;
typedef Arrangement_2::Halfedge_handle Halfedge_handle;


namespace barcode_templates {

  typedef int index;

  template<typename GrMat> void get_grades_of_generators(GrMat& M, std::vector<typename GrMat::Grade>& grades) {
    typedef typename GrMat::Grade Grade;
    double curr_x=std::numeric_limits<double>::min();
    double curr_y=std::numeric_limits<double>::min();
    for(index i = 0; i < M.num_rows; i++) {
      Grade& gr=M.row_grades[i];
      if(gr.at[0]==curr_x && gr.at[1]==curr_y) {
	continue;
      }
      grades.push_back(gr);
      curr_x=gr.at[0];
      curr_y=gr.at[1];
    }
  }

  template<typename GrMat> void get_grades_of_relations(GrMat& M, std::vector<typename GrMat::Grade>& grades) {
    typedef typename GrMat::Grade Grade;
    double curr_x=std::numeric_limits<double>::min();
    double curr_y=std::numeric_limits<double>::min();
    for(index i = 0; i < M.get_num_cols(); i++) {
      Grade& gr=M.grades[i];
      if(gr.at[0]==curr_x && gr.at[1]==curr_y) {
	continue;
      }
      grades.push_back(gr);
      curr_x=gr.at[0];
      curr_y=gr.at[1];
    }
  }

  template<typename Grade> bool are_incomparable(Grade& g1, Grade& g2) {
    return (g1.at[0] < g2.at[0] && g1.at[1] > g2.at[1]) || (g1.at[0] > g2.at[0] && g1.at[1] < g2.at[1]);
  }


  template<typename Grade> void get_joins(std::vector<Grade> grades, std::set<std::pair<typename Grade::Coordinate, typename Grade::Coordinate>>& joins) {

    index n = grades.size();
    for(index i=0;i<n;i++) {
      Grade& g1=grades[i];
      for(index j=i+1;j<n;j++) {
	Grade& g2=grades[j];
	if(are_incomparable(g1,g2)) {
	  //std::cout <<"Incomp: " << g1.first_val << " " << g1.second_val << " -- " << g2.first_val << " " << g2.second_val << std::endl;
	  joins.insert(std::make_pair(std::max(g1.at[0],g2.at[0]), std::max(g1.at[1],g2.at[1])));
	} else {
	  //std::cout <<"Not Incomp: " << g1.first_val << " " << g1.second_val << " -- " << g2.first_val << " " << g2.second_val << std::endl;
	}
      }
    }
  }


  template<typename Grade> void build_arrangement(Arrangement_2& arr,
						  std::set<std::pair<typename Grade::Coordinate, typename Grade::Coordinate>>& gen_joins,
						  std::set<std::pair<typename Grade::Coordinate, typename Grade::Coordinate>>& rel_joins) {
    

    typedef std::pair<typename Grade::Coordinate, typename Grade::Coordinate> Point;
    std::vector<X_monotone_curve_2> curves;
    for(Point p : gen_joins) {
      curves.push_back(Ray_2(Point_2(0,p.second),Point_2(1,p.second+p.first)));
    }
    for(Point p : rel_joins) {
      curves.push_back(Ray_2(Point_2(0,p.second),Point_2(1,p.second+p.first)));
    }
    curves.push_back(Line_2(Point_2(0,0),Point_2(0,1)));
    CGAL::insert(arr,curves.begin(),curves.end());
  }

  void add_joins_from_grades(const std::vector<matching::Grade2>& pts,
                             std::set<std::pair<double, double> >& joins) {
    std::vector<matching::Grade2> grades;
    std::set<matching::Grade2> seen;
    for (const matching::Grade2& g : pts) {
      if (seen.insert(g).second)
        grades.push_back(g);
    }
    const index n = static_cast<index>(grades.size());
    for (index i = 0; i < n; ++i) {
      for (index j = i + 1; j < n; ++j) {
        const matching::Grade2& g1 = grades[i];
        const matching::Grade2& g2 = grades[j];
        const bool incomparable =
            (g1.x < g2.x && g1.y > g2.y) || (g1.x > g2.x && g1.y < g2.y);
        if (!incomparable)
          continue;
        const double x = CGAL::to_double(g1.x > g2.x ? g1.x : g2.x);
        const double y = CGAL::to_double(g1.y > g2.y ? g1.y : g2.y);
        joins.insert(std::make_pair(x, y));
      }
    }
  }

  void add_joins_from_presentation(const matching::Presentation& P,
                                   std::set<std::pair<double, double> >& gen_joins,
                                   std::set<std::pair<double, double> >& rel_joins) {
    add_joins_from_grades(P.generators, gen_joins);
    add_joins_from_grades(P.relations, rel_joins);
  }

}


int main(int argc, char** argv) {

  if(argc<2) {
    std::cerr << "Usage: " << argv[0]
              << " file.scc|file.scc_k|file.sccsum|file.bif [file2] [--dim d] [--no-cgal] [--no-refine] [--no-aida]"
              << std::endl;
    std::exit(1);
  }

  std::vector<std::string> input_files;
  bool build_arrangement=true;
  bool refine=true;
  bool preprocess=true;
  int homology_dim = 0;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--no-cgal") {
      std::cout << "not building arrangement, only join lines" << std::endl;
      build_arrangement = false;
    } else if (arg == "--no-refine") {
      refine = false;
    } else if (arg == "--no-aida") {
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
  
  std::set<std::pair<double, double> > gen_joins;
  std::set<std::pair<double, double> > rel_joins;
  std::vector<matching::Module> modules;

  try {
  for (const std::string& input_filename : input_files) {
      matching::Module M = scc_hera::load_module(input_filename, homology_dim, preprocess);
      for (const matching::Presentation& P : M.summands)
          barcode_templates::add_joins_from_presentation(P, gen_joins, rel_joins);
      modules.push_back(std::move(M));
  }
  } catch (const std::exception& ex) {
      std::cerr << ex.what() << std::endl;
      std::exit(1);
  }

  // One file with exactly two presentations: treat them as the two modules
  // (old scc2020 two-block convention).
  if (modules.size() == 1 && modules[0].summands.size() == 2) {
      matching::Module A = matching::as_module(modules[0].summands[0]);
      matching::Module B = matching::as_module(modules[0].summands[1]);
      modules.clear();
      modules.push_back(std::move(A));
      modules.push_back(std::move(B));
  }

  std::cout << "In total, found " << gen_joins.size() << " generator joins and " << rel_joins.size() << " relation joins " << std::endl;
  std::cout << "Sum: " << gen_joins.size() + rel_joins.size() << std::endl;
  std::cout << "Loaded " << modules.size() << " module(s)";
  if (!modules.empty()) {
      std::cout << " with";
      for (std::size_t i = 0; i < modules.size(); ++i) {
          if (i) std::cout << " and";
          std::cout << " " << modules[i].summands.size() << " summand(s)";
      }
  }
  std::cout << std::endl;

  if(build_arrangement) {
      if (modules.size() != 2) {
          std::cerr << "Matching distance needs two modules (two files, or one file with two scc2020 blocks)" << std::endl;
          std::exit(1);
      }

      Arrangement_2 arr;
      matching::ArrangementStats stats;
      matching::build_arrangement_A(arr, modules[0], modules[1], refine, &stats);

      std::cout << "Arrangement has " << arr.number_of_vertices() << " vertices, "
                << arr.number_of_edges() << " edges, and "
                << arr.number_of_faces() << " faces" << std::endl;

      const matching::FT d = matching::maximize_F(arr, modules[0], modules[1], refine, &stats);
      std::cout << "Matching distance: ";
      if (matching::is_infinite_F(d))
          std::cout << "infinity" << std::endl;
      else
          std::cout << CGAL::to_double(d) << std::endl;
      std::cout << "A0 lines: " << stats.a0_lines << std::endl;
      std::cout << "Total lines: " << stats.lines_introduced << std::endl;
      std::cout << "Green-green deletions: " << stats.green_green_deletions << std::endl;
  }

  return 0;
}
