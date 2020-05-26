add_executable(testwrapper meemum/wrapper.test.cpp)
target_include_directories(testwrapper PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(testwrapper meemum gtest gtest_main)

#gtest_discover_tests(testwrapper WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/data/hgp)
add_test(NAME testwrapper COMMAND testwrapper 
         WORKING_DIRECTORY${PROJECT_SOURCE_DIR}/tests/testwrapper)
