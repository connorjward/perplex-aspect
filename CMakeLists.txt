# set minimum cmake version
cmake_minimum_required(VERSION 3.12)

# set the project name
project(PerplexAspect)

# enable Fortran compilation
enable_language(Fortran)

# enable testing
enable_testing()

add_subdirectory(src)
add_subdirectory(tests)

