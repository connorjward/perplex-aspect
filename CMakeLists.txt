# set minimum cmake version
cmake_minimum_required(VERSION 3.11)

# set the project name
project(PerpleXASPECT)

# required for tests to be found by ctest
# also declares BUILD_TESTING option
include(CTest)

# enable Fortran compilation
enable_language(Fortran)

set(CMAKE_CXX_STANDARD 11)

# enable position independent code (for shared libraries)
# TODO: make this target dependent?
# TODO: is this necessary?
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include_directories(include)

add_subdirectory(src)

