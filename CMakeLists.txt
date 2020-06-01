# set minimum cmake version
cmake_minimum_required(VERSION 3.11)

# set the project name
project(PerpleXASPECT)

# include modules
include(CTest)
include(FetchContent)
include(GoogleTest)

# enable Fortran compilation
enable_language(Fortran)

set(CMAKE_CXX_STANDARD 11)

# enable position independent code (for shared libraries)
# TODO: make this target dependent?
# TODO: is this necessary?
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(BUILD_TESTING)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG release-1.8.1
  )

  FetchContent_GetProperties(googletest)

  if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
  endif()
endif()

include_directories(include)
add_subdirectory(src)

