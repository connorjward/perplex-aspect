include(FetchContent)
include(GoogleTest)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.8.0
)

FetchContent_GetProperties(googletest)

if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
  add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif()

add_executable(meemum-test meemum/wrapper.test.cpp)

target_include_directories(
  meemum-test 
  PUBLIC ${PerpleXASPECT_SOURCE_DIR}/include
)

target_link_libraries(meemum-test meemum gtest_main)

# Copy test files to build directory
file(COPY ${PerpleXASPECT_SOURCE_DIR}/data/hgp/
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
     FILES_MATCHING PATTERN *.dat)

add_test(NAME meemum-test COMMAND meemum-test)
