include(FindGTest)
include(GoogleTest)

add_executable(khgp khgp/khgp.cpp)
target_include_directories(khgp PRIVATE ${PROJECT_SOURCE_DIR}/src ${GTEST_INCLUDE_DIRS})

target_link_libraries(khgp meemum_wrapper ${GTEST_BOTH_LIBRARIES})

gtest_discover_tests(khgp WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/klb1)
