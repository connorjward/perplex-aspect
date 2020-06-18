add_library(phaseinfo SHARED composition.cc)

target_link_libraries(phaseinfo perplexcpp)

ASPECT_SETUP_PLUGIN(phaseinfo)
