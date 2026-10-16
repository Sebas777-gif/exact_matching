#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include <mpp_utils/create_graded_matrices_from_scc2020.h>
#include <mpp_utils/Graded_matrix.h>
#include <scc/Scc.h>

#include <matching_distance.h>

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
              << " file.scc [file2.scc] [--no-cgal] [--no-refine]" << std::endl;
    std::exit(1);
  }

  std::vector<std::string> input_files;
  bool build_arrangement=true;
  bool refine=true;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--no-cgal") {
      std::cout << "not building arrangement, only join lines" << std::endl;
      build_arrangement = false;
    } else if (arg == "--no-refine") {
      refine = false;
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
  typedef GrMat::Grade Grade;

  std::vector<std::string> scc_strings;
  for (const std::string& input_filename : input_files) {
    split_file_into_scc2020_chunks(input_filename, scc_strings);
  }
  
  /*
  for(std::string s: scc_strings) {
      std::cout << "New string: \n\n";
      std::cout << "***" << s << "***" << std::endl;
  }
  */

  
  std::set<Point> gen_joins;
  std::set<Point> rel_joins;
  std::vector<matching::Presentation> modules;

  //scc::verbose=true;
  
  for(std::string s : scc_strings) {

      std::stringstream sstr(s);
      scc::Scc<Grade_cast> parser(sstr,&Grade_cast::get_instance()) ;
      
      std::vector<GrMat> matrices;
      mpp_utils::create_graded_matrices_from_scc2020(parser,matrices);
      
      GrMat& Mat = matrices[0];
      
      std::vector<Grade> gen_grades;
      barcode_templates::get_grades_of_generators<GrMat>(Mat,gen_grades);
      
      std::vector<Grade> rel_grades;
      barcode_templates::get_grades_of_relations<GrMat>(Mat,rel_grades);
      
      barcode_templates::get_joins(gen_grades,gen_joins);
      barcode_templates::get_joins(rel_grades,rel_joins);

      modules.push_back(barcode_templates::to_presentation(Mat));
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
      std::cout << "Matching distance: " << CGAL::to_double(d) << std::endl;
  }

  return 0;
}
