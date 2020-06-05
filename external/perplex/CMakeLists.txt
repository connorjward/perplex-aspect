add_library(
  perplex-obj 
  OBJECT
  BLASlib.f
  flib.f
  nlib.f
  olib.f
  resub.f
  rlib.f
  tlib.f
)

target_compile_options(perplex-obj PRIVATE -cpp)
target_compile_definitions(perplex-obj PRIVATE dgemv=x0001)
		       
