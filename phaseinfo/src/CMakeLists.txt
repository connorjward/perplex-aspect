add_library(phaseinfo SHARED 
            phase_composition.cc)

target_compile_options(phaseinfo PRIVATE -fext-numeric-literals)
target_link_libraries(phaseinfo perplex)

ASPECT_SETUP_PLUGIN(phaseinfo)
