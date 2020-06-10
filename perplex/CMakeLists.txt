cmake_minimum_required(VERSION 3.11)

project(PerpleX LANGUAGES CXX Fortran)

option(PERPLEX_BUILD_TESTING "docstring")

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
endif()

set(CMAKE_CXX_STANDARD 11) # TODO remove this

add_subdirectory(extern/perplex)
add_subdirectory(src)

# build the tests if this is the main CMake project 
if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR PERPLEX_BUILD_TESTING) AND BUILD_TESTING)
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

  add_subdirectory(test)
endif()
