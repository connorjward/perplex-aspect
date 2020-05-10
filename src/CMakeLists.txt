add_library(libmeemum meemum/interface.f meemum/wrapper.cpp)

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

