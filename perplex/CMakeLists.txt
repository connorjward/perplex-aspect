# set minimum cmake version
cmake_minimum_required(VERSION 3.11)

# set the project name
project(PerpleX LANGUAGES CXX Fortran)

# required for tests to be found by ctest
# also declares BUILD_TESTING option
include(CTest)

set(CMAKE_CXX_STANDARD 11)

# install GoogleTest
if(BUILD_TESTING)
  include(FetchContent)
  include(GoogleTest)

  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        release-1.8.0
  )

  FetchContent_GetProperties(googletest)

  if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
  endif()
endif()

add_subdirectory(extern/perplex)
add_subdirectory(src)

