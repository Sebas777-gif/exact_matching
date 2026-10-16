#include <iostream>

#include <mpp_utils/create_graded_matrices_from_scc2020.h>
#include <mpp_utils/Graded_matrix.h>
#include <scc/Scc.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Arr_linear_traits_2<Kernel>             Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Segment_2                           Segment_2;
typedef Traits_2::Ray_2                               Ray_2;
typedef Traits_2::Line_2                              Line_2;
typedef Traits_2::X_monotone_curve_2                  X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef Arrangement_2::Vertex_handle                  Vertex_handle;
typedef Arrangement_2::Halfedge_handle                Halfedge_handle;


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
    std::cerr << "Argument!" << std::endl;
    std::exit(1);
  }

  std::string input_filename=argv[1];

  bool build_arrangement=true;

  if(argc>=3 && std::string(argv[2])=="--no-cgal") {
      std::cout << "not building arrangment, only join lines" << std::endl;
      build_arrangement=false;
  }
  
  typedef mpp_utils::Graded_matrix<> GrMat;
  typedef GrMat::Coordinate Coordinate;
  typedef std::pair<Coordinate,Coordinate> Point;
  
  typedef GrMat::Coordinate_traits Grade_cast;
  typedef GrMat::Grade Grade;
  
  std::vector<std::vector<GrMat> > scc_instances;

  std::vector<std::string> scc_strings;

  split_file_into_scc2020_chunks(input_filename, scc_strings);
  
  /*
  for(std::string s: scc_strings) {
      std::cout << "New string: \n\n";
      std::cout << "***" << s << "***" << std::endl;
  }
  */

  
  std::set<Point> gen_joins;
  std::set<Point> rel_joins;

  //scc::verbose=true;
  
  for(std::string s : scc_strings) {

      std::stringstream sstr(s);
      scc::Scc<Grade_cast> parser(sstr,&Grade_cast::get_instance()) ;
      
      std::vector<GrMat> matrices;
      mpp_utils::create_graded_matrices_from_scc2020(parser,matrices);
      
      GrMat& M = matrices[0];
      
      
      //std::cout << "Found scc2020 instance" << std::endl;
      //std::cout << M.get_num_cols() << " " << M.num_rows << std::endl;
      
      std::vector<Grade> gen_grades;
      barcode_templates::get_grades_of_generators<GrMat>(M,gen_grades);
      
      //std::cout << "Found " << gen_grades.size() << " different generator grades" << std::endl;
      
      std::vector<Grade> rel_grades;
      barcode_templates::get_grades_of_relations<GrMat>(M,rel_grades);
      
      //std::cout << "Found " << rel_grades.size() << " different relation grades" << std::endl;
      
      
      
      barcode_templates::get_joins(gen_grades,gen_joins);
      
      //std::cout << "Found " << gen_joins.size() << " joins for generators (maximal is " << gen_grades.size()*(gen_grades.size()-1)/2 <<")" << std::endl;
      
      barcode_templates::get_joins(rel_grades,rel_joins);
      
      //std::cout << "Found " << rel_joins.size() << " joins for relations (maximal is " << rel_grades.size()*(rel_grades.size()-1)/2 <<")" << std::endl;
      
  }

  std::cout << "In total, found " << gen_joins.size() << " generator joins and " << rel_joins.size() << " relation joins " << std::endl;
  std::cout << "Sum: " << gen_joins.size() + rel_joins.size() << std::endl;

  if(build_arrangement) {
  
      Arrangement_2 arr;
      
      barcode_templates::build_arrangement<Grade>(arr,gen_joins,rel_joins);
      

      std::cout << "Arrangement has " << arr.number_of_vertices() << " vertices, " << arr.number_of_edges() << " edges, and " << arr.number_of_faces() << " faces" << std::endl;
  }
  /*
  for(Grade& g : gen_grades) {
    std::cout << g.first_val << " " << g.second_val << std::endl;
  }
  */

  return 0;
}
