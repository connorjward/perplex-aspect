add_executable(perplex-test 
               c_interface.cc state.cc)

target_include_directories(perplex-test PRIVATE ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(perplex-test perplex gtest gtest_main)

gtest_discover_tests(perplex-test
		     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test-simple)

# copy test files to build directory
file(COPY ${CMAKE_SOURCE_DIR}/data/perplex-dat-files/test-simple
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
     FILES_MATCHING PATTERN *.dat)
