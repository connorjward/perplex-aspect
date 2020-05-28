add_library(meemum wrapper.cc ftoc.f)

# include Perple_X header files
target_include_directories(meemum 
                           PUBLIC ${PerpleXASPECT_SOURCE_DIR}/include
                           PRIVATE $ENV{PERPLEX_DIR})

# link with precompiled Perple_X object code
target_link_libraries(meemum 
                      PRIVATE $ENV{PERPLEX_DIR}/resub.o 
                              $ENV{PERPLEX_DIR}/flib.o 
                              $ENV{PERPLEX_DIR}/nlib.o 
                              $ENV{PERPLEX_DIR}/olib.o 
                              $ENV{PERPLEX_DIR}/rlib.o 
                              $ENV{PERPLEX_DIR}/tlib.o)
