# Copyright (C) 2011 - 2016 by the authors of the ASPECT code.
# Copyright (C) 2020 by Connor Ward.

cmake_minimum_required(VERSION 3.11)

project(PerpleX-ASPECT CXX)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
endif()

add_subdirectory(src)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  # copy test files
  file(COPY ${PROJECT_SOURCE_DIR}/data/test
       DESTINATION ${PROJECT_BINARY_DIR}/data)

  add_test(NAME CheckRuns
           COMMAND ${PROJECT_BINARY_DIR}/src/aspect test.prm
	   WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/data/test/simple)
endif()
