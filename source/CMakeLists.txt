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

# add_library(perplexaspect SHARED perplex_melt_global.cc perplex_particle.cc)
add_library(perplexaspect SHARED perplex_melt_simple.cc)
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
