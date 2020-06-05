# library
add_library(
  perplex
  interface.f
  solver.cc 
)

target_include_directories(
  perplex 
  PRIVATE ${PerpleXASPECT_SOURCE_DIR}/external/perplex
)

target_link_libraries(
  perplex
  PRIVATE perplex-obj
)
		       
# add_custom_command(
#   TARGET perplex
#   PRE_BUILD
#   COMMAND sed -i 's/dgemv/x0001/g' ${perplex-tmp_SOURCE_DIR}/BLASlib.f
#   COMMAND sed -i 's/dgemv/x0001/g' ${perplex-tmp_SOURCE_DIR}/nlib.f
#   COMMENT "Replace all occurences of dgemv with x0001 to prevent linking conflicts with BLAS."
# )

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
