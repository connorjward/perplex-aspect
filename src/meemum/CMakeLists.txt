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

add_executable(meemum-test wrapper.test.cc)

target_include_directories(
  meemum-test 
  PUBLIC ${PerpleXASPECT_SOURCE_DIR}/include
)

target_link_libraries(meemum-test meemum gtest_main)

# Copy test files to build directory
file(
  COPY ${PerpleXASPECT_SOURCE_DIR}/data/hgp
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME meemum-test COMMAND meemum-test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hgp
)
