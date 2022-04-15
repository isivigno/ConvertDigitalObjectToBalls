# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12)
project(convertDigitalUoB_5 CXX)

set(CMAKE_BUILD_TYPE Release)	

set(CMAKE_CXX_FLAGS "-O3")

#include for local directory
############################
file(
        GLOB
        source_files
        src/*.h src/*.ih src/convertDigitalUoB.cpp
)


# ##################################################

add_executable(convertDigitalUoB  ${source_files})	

			       
################## External libraries #################

## DGtal
########
FIND_PACKAGE(DGtal REQUIRED)

INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

## OpenMP
#########
find_package(OpenMP REQUIRED)


## CGAL and its components
##########################
find_package( CGAL QUIET COMPONENTS Core Qt5)

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

## Boost and its components
###########################
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()


###########################################################

target_link_libraries(convertDigitalUoB  PRIVATE ${DGTAL_LIBRARIES} ${UNIONOFDISKS} OpenMP::OpenMP_CXX CGAL::CGAL_Core ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} CGAL::CGAL_Qt5 )


