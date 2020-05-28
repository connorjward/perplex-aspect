add_executable(meemum-test meemum/wrapper.test.cpp)
target_include_directories(meemum-test 
			   PUBLIC ${PerpleXASPECT_SOURCE_DIR}/include
                           PRIVATE ${PerpleXASPECT_SOURCE_DIR}/src)
target_link_libraries(meemum-test meemum gtest_main)

# Copy test files to build directory
file(COPY ${PerpleXASPECT_SOURCE_DIR}/data/hgp/
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
     FILES_MATCHING PATTERN *.dat)

add_test(NAME meemum-test COMMAND meemum-test)
