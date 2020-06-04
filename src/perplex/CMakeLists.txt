#temp
include(FetchContent)
FetchContent_Declare(
  perplex-tmp
  URL https://petrol.natur.cuni.cz/~ondro/perplex-sources-stable.zip
)

FetchContent_GetProperties(perplex-tmp)

if(NOT perplex-tmp_POPULATED)
  FetchContent_Populate(perplex-tmp)
endif()

# library
add_library(
  perplex
  solver.cc 
  interface.f
  ${perplex-tmp_SOURCE_DIR}/BLASlib.f
  ${perplex-tmp_SOURCE_DIR}/flib.f
  ${perplex-tmp_SOURCE_DIR}/nlib.f
  ${perplex-tmp_SOURCE_DIR}/olib.f
  ${perplex-tmp_SOURCE_DIR}/resub.f
  ${perplex-tmp_SOURCE_DIR}/rlib.f
  ${perplex-tmp_SOURCE_DIR}/tlib.f
)

target_include_directories(
  perplex 
  PUBLIC ${PerpleXASPECT_SOURCE_DIR}/include
  PRIVATE ${perplex-tmp_SOURCE_DIR}
)
		       
add_custom_command(
  TARGET perplex
  PRE_BUILD
  COMMAND sed -i 's/dgemv/x0001/g' ${perplex-tmp_SOURCE_DIR}/BLASlib.f
  COMMAND sed -i 's/dgemv/x0001/g' ${perplex-tmp_SOURCE_DIR}/nlib.f
  COMMENT "Replace all occurences of dgemv with x0001 to prevent linking conflicts with BLAS."
)

# unit tests
add_executable(
  test-perplex
  interface.test.cc
  solver.test.cc
)

target_link_libraries(
  test-perplex
  perplex
  gtest
  gtest_main
)

# copy test files to build directory
file(
  COPY ${PerpleXASPECT_SOURCE_DIR}/data/perplex/tests
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)

# automatically add tests
gtest_discover_tests(
  test-perplex
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests
)
