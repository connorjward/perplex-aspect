# set minimum cmake version
cmake_minimum_required(VERSION 3.6)

# set the project name
project(PerpleXASPECT)

# Enable Fortran compilation
enable_language(Fortran)

# Enable position independent code (for shared libraries)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CTest)
include(FetchContent)
include(GoogleTest)

if(BUILD_TESTING)
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

add_subdirectory(src)

