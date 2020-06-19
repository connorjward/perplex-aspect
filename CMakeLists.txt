# Copyright (C) 2011 - 2016 by the authors of the ASPECT code.
# Copyright (C) 2020 by Connor Ward.

cmake_minimum_required(VERSION 3.11)

project(aspect-phaseinfo CXX)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
endif()

find_package(Aspect 2.1.0 QUIET HINTS ${Aspect_DIR} ../ ../../ $ENV{ASPECT_DIR})

if(NOT Aspect_FOUND)
  message(FATAL_ERROR "\n"
	"Could not find a valid ASPECT build/installation directory. "
	"Please specify the directory where you are building ASPECT by passing\n"
	"   -D Aspect_DIR=<path to ASPECT>\n"
	"to cmake or by setting the environment variable ASPECT_DIR in your shell "
	"before calling cmake. See the section 'How to write a plugin' in the "
        "manual for more information.")
endif()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

include(FetchContent)
FetchContent_Declare(
  perplexcpp
  GIT_REPOSITORY https://github.com/cward97/perplex-cpp.git
  GIT_TAG        origin/master
)

FetchContent_GetProperties(perplexcpp)
if(NOT perplexcpp_POPULATED)
  FetchContent_Populate(perplexcpp)
  add_subdirectory(${perplexcpp_SOURCE_DIR} ${perplexcpp_BINARY_DIR})
endif()

add_subdirectory(src)

# Build the tests if this is the main CMake project 
# or if PHASEINFO_BUILD_TESTING is enabled.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  # copy test files
  file(COPY ${PROJECT_SOURCE_DIR}/data/test
       DESTINATION ${PROJECT_BINARY_DIR}/data)

  add_test(NAME CheckRuns
           COMMAND ${PROJECT_BINARY_DIR}/src/aspect test.prm
	   WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/data/test/simple)
endif()
