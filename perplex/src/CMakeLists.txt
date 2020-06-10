add_library(perplex SHARED
	    c_interface.f meemum_solver.cc)

target_include_directories(perplex 
			   PUBLIC ${PROJECT_SOURCE_DIR}/include
			   PRIVATE ${perplex-extern_SOURCE_DIR})

target_link_libraries(perplex PRIVATE perplex-extern)

target_compile_options(perplex
		       PRIVATE 
		       -ffixed-line-length-132) # increase line length
