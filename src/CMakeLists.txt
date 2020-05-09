add_library(libmeemum libmeemum.f)

target_include_directories(libmeemum PRIVATE $ENV{PERPLEX_DIR})

set(perplex_libs resub.o flib.o nlib.o olib.o rlib.o tlib.o)
list(TRANSFORM perplex_libs PREPEND $ENV{PERPLEX_DIR}/)
target_link_libraries(libmeemum PRIVATE ${perplex_libs})

add_executable(decompr-melt decompr_melt.cc meemum_wrapper.cc)

target_link_libraries(decompr-melt libmeemum)

