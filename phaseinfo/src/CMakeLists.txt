add_library(phaseinfo SHARED composition.cc)

target_link_libraries(phaseinfo perplex)

ASPECT_SETUP_PLUGIN(phaseinfo)
