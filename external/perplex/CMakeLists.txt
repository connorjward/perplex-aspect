add_library(
  perplex-extern
  OBJECT
  BLASlib.f
  flib.f
  nlib.f
  olib.f
  resub.f
  rlib.f
  tlib.f
)

set_target_properties(perplex-extern PROPERTIES POSITION_INDEPENDENT_CODE True)

target_compile_options(
  perplex-extern 
  PRIVATE 
  -cpp			  # enable preprocessing 
  -ffixed-line-length-132 # increase line length so prefix doesn't break compilation
)

# add a prefix to conflicting function names
target_compile_definitions(
  perplex-extern
  PRIVATE 
  dgemv=perplex_dgemv # conflicts with BLAS dgemv
)
		       
set(perplex-extern_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} 
    CACHE PATH "Location of perplex-extern source files")
