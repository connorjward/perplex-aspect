# Copyright (C) 2015 - 2020 by the authors of the ASPECT code.
# Copyright (C) 2020 Connor Ward
#
# This file is part of PerpleX-ASPECT.
# PerpleX-ASPECT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PerpleX-ASPECT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PerpleX-ASPECT. If not, see <https://www.gnu.org/licenses/>.


# Find ASPECT.
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


# Build the main library.
add_library(
  perplexaspect 
  SHARED 
  initial_composition/perplex_composition.cc
  material_model/perplex_melt.cc 
  particle/property/perplex_composition.cc
  postprocess/perplex_cache_statistics.cc
)

target_link_libraries(perplexaspect perplexcpp)

ASPECT_SETUP_PLUGIN(perplexaspect)


# To allow the cookbooks to be straightforwardly run, the ASPECT executable 
# (linked to the current directory via ASPECT_SETUP_PLUGIN) and phaseinfo 
# shared library are symlinked into the cookbooks directory.
set(COOKBOOKS_DIR ${PROJECT_SOURCE_DIR}/cookbooks)

add_custom_command(
  OUTPUT ${COOKBOOKS_DIR}/aspect ${COOKBOOKS_DIR}/libperplexaspect.so
  COMMAND ${CMAKE_COMMAND} -E create_symlink 
	  ${CMAKE_CURRENT_BINARY_DIR}/aspect ${COOKBOOKS_DIR}/aspect
  COMMAND ${CMAKE_COMMAND} -E create_symlink 
	  ${CMAKE_CURRENT_BINARY_DIR}/libperplexaspect.so ${COOKBOOKS_DIR}/libperplexaspect.so)

add_custom_target(setup_cookbooks
  DEPENDS ${COOKBOOKS_DIR}/aspect ${COOKBOOKS_DIR}/libperplexaspect.so)

add_dependencies(setup_cookbooks perplexaspect)


# Do the same for the tests.
set(TESTS_DIR ${PROJECT_SOURCE_DIR}/tests)

add_custom_command(
  OUTPUT ${TESTS_DIR}/aspect ${TESTS_DIR}/libperplexaspect.so
  COMMAND ${CMAKE_COMMAND} -E create_symlink 
	  ${CMAKE_CURRENT_BINARY_DIR}/aspect ${TESTS_DIR}/aspect
  COMMAND ${CMAKE_COMMAND} -E create_symlink 
	  ${CMAKE_CURRENT_BINARY_DIR}/libperplexaspect.so ${TESTS_DIR}/libperplexaspect.so)

add_custom_target(setup_tests
  DEPENDS ${TESTS_DIR}/aspect ${TESTS_DIR}/libperplexaspect.so)

add_dependencies(setup_tests perplexaspect)
