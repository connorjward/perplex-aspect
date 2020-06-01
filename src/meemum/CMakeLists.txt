
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

target_include_directories(
  meemum 
  PUBLIC ${PerpleXASPECT_SOURCE_DIR}/include
  PRIVATE ${perplex_SOURCE_DIR}
)
		       
add_custom_command(
  TARGET meemum
  PRE_BUILD
  COMMAND sed -i 's/dgemv/x0001/g' ${perplex_SOURCE_DIR}/BLASlib.f
  COMMAND sed -i 's/dgemv/x0001/g' ${perplex_SOURCE_DIR}/nlib.f
  COMMENT "Replace all occurences of dgemv with x0001 to prevent linking conflicts with BLAS."
)

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
