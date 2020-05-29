# set minimum cmake version
cmake_minimum_required(VERSION 3.6)

# set the project name
project(PerpleXASPECT)

# Enable Fortran compilation
enable_language(Fortran)

# Enable position independent code (for shared libraries)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


add_subdirectory(src)

option(PACKAGE_TESTS "Build the tests" ON)
if(PACKAGE_TESTS)
  enable_testing()

  add_subdirectory(tests)
endif()


