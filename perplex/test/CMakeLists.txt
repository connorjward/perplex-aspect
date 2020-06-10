add_executable(perplex_test c_interface.cc meemum_solver.cc)

target_include_directories(perplex_test PRIVATE ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(perplex_test perplex gtest gtest_main)

gtest_discover_tests(perplex_test
		     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test-fast)

# copy test files to build directory
file(COPY ${CMAKE_SOURCE_DIR}/data/perplex-dat-files/test-fast
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
     FILES_MATCHING PATTERN *.dat)
