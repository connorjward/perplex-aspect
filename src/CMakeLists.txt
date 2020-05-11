add_library(libmeemum meemum/calcs.f meemum/props.f meemum/wrapper.cpp)

target_include_directories(libmeemum 
			   PUBLIC ../include
			   PRIVATE $ENV{PERPLEX_DIR})

set(perplex_libs resub.o flib.o nlib.o olib.o rlib.o tlib.o)
list(TRANSFORM perplex_libs PREPEND $ENV{PERPLEX_DIR}/)
target_link_libraries(libmeemum PRIVATE ${perplex_libs})

add_executable(decompr-melt decompr_melt.cpp)

target_include_directories(decompr-melt 
			   PRIVATE ../include)

target_link_libraries(decompr-melt libmeemum)

# testing
add_executable(misc misc.cpp)

target_link_libraries(misc libmeemum)

include(FindGTest)
include(GoogleTest)

add_executable(wrapper-test meemum/wrapper.test.cpp)
target_include_directories(wrapper-test PRIVATE ${PROJECT_SOURCE_DIR}/include ${GTEST_INCLUDE_DIRS})

target_link_libraries(wrapper-test libmeemum ${GTEST_BOTH_LIBRARIES})

gtest_discover_tests(wrapper-test WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/data/hgp)
