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

include_directories(include)

add_subdirectory(external/perplex)

add_subdirectory(src/aspect)
add_subdirectory(src/perplex)

