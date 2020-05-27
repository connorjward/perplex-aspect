add_library(meemum meemum/meemum.f meemum/props.f)

# include Perple_X header files
target_include_directories(meemum PRIVATE $ENV{PERPLEX_DIR})

# link with precompiled Perple_X object code
target_link_libraries(meemum 
                      PRIVATE $ENV{PERPLEX_DIR}/resub.o 
                              $ENV{PERPLEX_DIR}/flib.o 
                              $ENV{PERPLEX_DIR}/nlib.o 
                              $ENV{PERPLEX_DIR}/olib.o 
                              $ENV{PERPLEX_DIR}/rlib.o 
                              $ENV{PERPLEX_DIR}/tlib.o)

######

add_executable(decompr-melt decompr_melt.cpp)

target_link_libraries(decompr-melt meemum)

add_subdirectory(plugin)
