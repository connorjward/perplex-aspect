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

target_compile_options(
  perplex
  PRIVATE 
  -ffixed-line-length-132 # increase line length 
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
  COPY ${PerpleXASPECT_SOURCE_DIR}/data/perplex-dat-files/test1
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  FILES_MATCHING PATTERN *.dat
)

# automatically add tests
gtest_discover_tests(
  test-perplex
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test1
)
