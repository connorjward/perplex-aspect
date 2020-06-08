# library
add_library(
  perplex
  SHARED
  interface.f
  solver.cc 
)

target_include_directories(
  perplex 
  PRIVATE ${perplex-extern_SOURCE_DIR}
)

target_link_libraries(
  perplex
  PRIVATE perplex-extern
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
