#include(FindBLAS)
#include(FindLAPACK)
#
#if(NOT BLAS_FOUND)
#  # some error checking...
#  message("BLAS library not found")
#else()
#  message("Using BLAS libraries found at: ${BLAS_LIBRARIES}")
#endif()
#
#if(NOT LAPACK_FOUND)
#  # some error checking...
#  message("LAPACK library not found")
#else()
#  message("Using LAPACK libraries found at: ${LAPACK_LIBRARIES}")
#endif()

FetchContent_Declare(
  perplex
  URL https://petrol.natur.cuni.cz/~ondro/perplex-sources-stable.zip
)

FetchContent_GetProperties(perplex)

if(NOT perplex_POPULATED)
  FetchContent_Populate(perplex)
endif()

add_library(
  meemum 
  wrapper.cc 
  ftoc.f
  ${perplex_SOURCE_DIR}/BLASlib.f
  ${perplex_SOURCE_DIR}/flib.f
  ${perplex_SOURCE_DIR}/nlib.f
  ${perplex_SOURCE_DIR}/olib.f
  ${perplex_SOURCE_DIR}/resub.f
  ${perplex_SOURCE_DIR}/rlib.f
  ${perplex_SOURCE_DIR}/tlib.f
)

# include Perple_X header files
target_include_directories(meemum 
                           PUBLIC ${PerpleXASPECT_SOURCE_DIR}/include
			   PRIVATE ${perplex_SOURCE_DIR})

#		 target_link_libraries(meemum ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# link with precompiled Perple_X object code
#target_link_directories(meemum ${perplex_SOURCE_DIR})
#target_link_libraries(meemum 
#                      PRIVATE $ENV{PERPLEX_DIR}/resub.o 
#                              $ENV{PERPLEX_DIR}/flib.o 
#                              $ENV{PERPLEX_DIR}/nlib.o 
#                              $ENV{PERPLEX_DIR}/olib.o 
#                              $ENV{PERPLEX_DIR}/rlib.o 
#                              $ENV{PERPLEX_DIR}/tlib.o)

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
