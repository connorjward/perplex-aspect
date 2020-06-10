add_library(phaseinfo SHARED phase_composition.cc)

target_link_libraries(phaseinfo perplex)

ASPECT_SETUP_PLUGIN(phaseinfo)
