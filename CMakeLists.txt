# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( cavity_volumes_pbc )

cmake_minimum_required(VERSION 2.8.11)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  
endif()

# include helper file
include( ${CGAL_USE_FILE} )

# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()  
endif()
include_directories( ${Boost_INCLUDE_DIRS} )

# include for local directory

# include for local package

# Creating entries for target: cavity_volumes_pbc
# ############################
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function -Wno-long-long")

add_executable( cavity_volumes_pbc  cavity_volumes_pbc.cpp config_file.cpp skin_surface.cpp sastry_subsimplex_volume.c )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS cavity_volumes_pbc )

# Link the executable to CGAL and third-party libraries
target_link_libraries(cavity_volumes_pbc   ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
