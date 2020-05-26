add_executable(testwrapper meemum/wrapper.test.cpp)
target_include_directories(testwrapper PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(testwrapper meemum gtest gtest_main)

# Copy test files to build directory
file(COPY ${PerplexAspect_SOURCE_DIR}/data/hgp/
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
     FILES_MATCHING PATTERN *.dat)

add_test(NAME testwrapper COMMAND testwrapper)
