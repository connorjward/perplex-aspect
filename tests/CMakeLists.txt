include(FindGTest)
include(GoogleTest)

add_executable(testwrapper meemum/wrapper.test.cpp)
target_include_directories(testwrapper PRIVATE ${PROJECT_SOURCE_DIR}/src ${GTEST_INCLUDE_DIRS})
target_link_libraries(testwrapper meemum ${GTEST_BOTH_LIBRARIES})

gtest_discover_tests(testwrapper WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/data/hgp)
