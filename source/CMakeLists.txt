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

add_library(perplexparticle SHARED perplex_particle.cc)
target_link_libraries(perplexparticle perplexcpp)
ASPECT_SETUP_PLUGIN(perplexparticle)

# To allow the cookbooks to be straightforwardly run, the ASPECT executable 
# (linked to the current directory via ASPECT_SETUP_PLUGIN) and phaseinfo 
# shared library are symlinked into the cookbooks directory.
set(COOKBOOKS_DIR ${PROJECT_SOURCE_DIR}/cookbooks)

add_custom_command(
  OUTPUT ${COOKBOOKS_DIR}/aspect ${COOKBOOKS_DIR}/libperplexparticle.so
  COMMAND ${CMAKE_COMMAND} -E create_symlink 
	  ${CMAKE_CURRENT_BINARY_DIR}/aspect ${COOKBOOKS_DIR}/aspect
  COMMAND ${CMAKE_COMMAND} -E create_symlink 
	  ${CMAKE_CURRENT_BINARY_DIR}/libperplexparticle.so ${COOKBOOKS_DIR}/libperplexparticle.so)

add_custom_target(setup_cookbooks
  DEPENDS ${COOKBOOKS_DIR}/aspect ${COOKBOOKS_DIR}/libperplexparticle.so)

add_dependencies(setup_cookbooks perplexparticle)
