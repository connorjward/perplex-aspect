add_library(meemum calcs.f props.f wrapper.cpp)

# include Perple_X header files
target_include_directories(meemum PRIVATE $ENV{PERPLEX_DIR})

# link with precompiled Perple_X object code
set(perplex_libs resub.o flib.o nlib.o olib.o rlib.o tlib.o)
list(TRANSFORM perplex_libs PREPEND $ENV{PERPLEX_DIR}/)

target_link_libraries(libmeemum PRIVATE ${perplex_libs})
