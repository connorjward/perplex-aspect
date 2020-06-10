add_library(phaseinfo SHARED 
            phase_composition.cc)

find_package(MPI REQUIRED)
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
target_link_libraries(phaseinfo MPI::MPI_CXX perplex)

ASPECT_SETUP_PLUGIN(phaseinfo)
