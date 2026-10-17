#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <mpp_utils/create_graded_matrices_from_scc2020.h>
#include <mpp_utils/Graded_matrix.h>
#include <scc/Scc.h>

#include <matching_distance.h>
#include <bif_io.hpp>

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

  typedef mpp_utils::index index;

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

  template<typename GrMat>
  matching::Presentation to_presentation(const GrMat& M) {
    matching::Presentation P;
    P.generators.resize(M.num_rows);
    for (index i = 0; i < M.num_rows; ++i) {
      P.generators[i].x = matching::FT(M.row_grades[i].at[0]);
      P.generators[i].y = matching::FT(M.row_grades[i].at[1]);
    }
    const index n_rel = M.get_num_cols();
    P.relations.resize(n_rel);
    P.columns.resize(n_rel);
    for (index i = 0; i < n_rel; ++i) {
      P.relations[i].x = matching::FT(M.grades[i].at[0]);
      P.relations[i].y = matching::FT(M.grades[i].at[1]);
      std::vector<index> col;
      M.get_col(i, col);
      P.columns[i].assign(col.begin(), col.end());
    }
    return P;
  }

}


void split_file_into_scc2020_chunks(std::string input_filename,std::vector<std::string>& scc_instances) {

    std::ifstream ifstr(input_filename);

    std::string next, curr_scc;

    
    while(ifstr.good()) {

	std::getline(ifstr,next);

	if(next=="") {
	    continue;
	}

	if(next=="scc2020") {
	    if(curr_scc!="") {
		scc_instances.push_back(curr_scc);
		curr_scc="";
	    }
	}
	curr_scc.append(next);
	curr_scc.append("\n");
    }
    if(curr_scc!="") {
	scc_instances.push_back(curr_scc);
    }
    
}


int main(int argc, char** argv) {

  if(argc<2) {
    std::cerr << "Usage: " << argv[0]
              << " file.scc|file.bif [file2.scc|file2.bif] [--dim d] [--no-cgal] [--no-refine]"
              << std::endl;
    std::exit(1);
  }

  std::vector<std::string> input_files;
  bool build_arrangement=true;
  bool refine=true;
  int homology_dim = 0;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--no-cgal") {
      std::cout << "not building arrangement, only join lines" << std::endl;
      build_arrangement = false;
    } else if (arg == "--no-refine") {
      refine = false;
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
  
  typedef mpp_utils::Graded_matrix<> GrMat;
  typedef GrMat::Coordinate Coordinate;
  typedef std::pair<Coordinate,Coordinate> Point;
  
  typedef GrMat::Coordinate_traits Grade_cast;

  std::set<Point> gen_joins;
  std::set<Point> rel_joins;
  std::vector<matching::Presentation> modules;

  auto take_presentation = [&](matching::Presentation P) {
      barcode_templates::add_joins_from_presentation(P, gen_joins, rel_joins);
      modules.push_back(std::move(P));
  };

  //scc::verbose=true;

  try {
  for (const std::string& input_filename : input_files) {
      if (scc_hera::is_bif_path(input_filename)) {
          take_presentation(scc_hera::load_bif_presentation(input_filename, homology_dim));
          continue;
      }

      std::vector<std::string> scc_strings;
      split_file_into_scc2020_chunks(input_filename, scc_strings);
      for (std::string s : scc_strings) {
          std::stringstream sstr(s);
          scc::Scc<Grade_cast> parser(sstr, &Grade_cast::get_instance());
          std::vector<GrMat> matrices;
          mpp_utils::create_graded_matrices_from_scc2020(parser, matrices);
          take_presentation(barcode_templates::to_presentation(matrices[0]));
      }
  }
  } catch (const std::exception& ex) {
      std::cerr << ex.what() << std::endl;
      std::exit(1);
  }

  std::cout << "In total, found " << gen_joins.size() << " generator joins and " << rel_joins.size() << " relation joins " << std::endl;
  std::cout << "Sum: " << gen_joins.size() + rel_joins.size() << std::endl;
  std::cout << "Loaded " << modules.size() << " presentation(s)" << std::endl;

  if(build_arrangement) {
      if (modules.size() < 2) {
          std::cerr << "Matching distance needs two presentations (two scc2020 blocks or two files)" << std::endl;
          std::exit(1);
      }

      Arrangement_2 arr;
      matching::build_arrangement_A(arr, modules[0], modules[1], refine);

      std::cout << "Arrangement has " << arr.number_of_vertices() << " vertices, "
                << arr.number_of_edges() << " edges, and "
                << arr.number_of_faces() << " faces" << std::endl;

      const matching::FT d = matching::maximize_F(arr, modules[0], modules[1], refine);
      std::cout << "Matching distance: ";
      if (matching::is_infinite_F(d))
          std::cout << "infinity" << std::endl;
      else
          std::cout << CGAL::to_double(d) << std::endl;
  }

  return 0;
}
