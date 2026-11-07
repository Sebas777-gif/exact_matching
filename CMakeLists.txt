cmake_minimum_required(VERSION 3.10)

project( Code_ )

set(CMAKE_CXX_STANDARD 17)



if (NOT CMAKE_VERBOSE_MAKEFILE)
   set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if (NOT PHAT_INCLUDE_DIR)
    set(PHAT_INCLUDE_DIR "~/phat/include")
endif()

if (NOT MPP_UTILS_INCLUDE_DIR)
    set(MPP_UTILS_INCLUDE_DIR "~/mpfree/mpp_utils_mod/include")
endif()

if (NOT SCC_INCLUDE_DIR)
    set(SCC_INCLUDE_DIR "~/chain_complex_format/include")
endif()

if (NOT CGAL_DIR)
   set(CGAL_DIR "~/CGAL-4.14")
endif()

include_directories(${PHAT_INCLUDE_DIR})
include_directories(${MPP_UTILS_INCLUDE_DIR})
include_directories(${SCC_INCLUDE_DIR})
include_directories(include)

#FIND_PACKAGE(OpenMP)

#find_package(Boost OPTIONAL_COMPONENTS timer)

find_package(CGAL QUIET COMPONENTS Core )

if(Boost_TIMER_FOUND)
    add_compile_options(-DTIMERS=1)
    link_libraries(Boost::timer)
endif()

if(CGAL_FOUND)
	add_executable(play "play.cpp")
	target_link_libraries(play ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
else()
        message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()
