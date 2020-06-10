cmake_minimum_required(VERSION 2.8)
project(PerpleX-ASPECT CXX Fortran)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
endif()

add_subdirectory(perplex)
add_subdirectory(phaseinfo)
