add_library(perplex SHARED
	    state.cc utils.cc c_interface.f)

target_include_directories(perplex 
			   PUBLIC ${PROJECT_SOURCE_DIR}/include
			   PRIVATE ${perplex-extern_SOURCE_DIR})

target_link_libraries(perplex PRIVATE perplex-extern)

target_compile_options(perplex
		       PRIVATE 
		       -ffixed-line-length-132) # increase line length
