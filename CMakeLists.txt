# set minimum cmake version
cmake_minimum_required(VERSION 3.12)
enable_language(Fortran)

# set the project name
project(PerplexAspect)

add_subdirectory(src)

enable_testing()
add_subdirectory(tests)
