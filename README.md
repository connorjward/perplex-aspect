# PerpleX-ASPECT

This repository contains two libraries: 

- `perplex` : A C++ wrapper for the thermodynamic code [Perple_X](perplex.ethz.ch) (written in Fortran 77).

- `phaseinfo` : A plugin for the geodynamical code [ASPECT](aspect.geodynamics.org) that can analyse phase (mineral) amounts and composition.

## Prerequisites

- [ASPECT](github.com/geodynamics/aspect)

- CMake version `3.11+`

- A Fortran 2003 compatible Fortran compiler (e.g. `gfortran`)

- A C++ 2011 compatible C++ compiler (e.g. `gcc`)

- A MPI library (e.g. `openmpi`)

## Installation instructions

An out of source build is recommended:

	mkdir build
	cd build
	cmake -DAspect_DIR=/path/to/aspect/ ..
	make -j<N>
	
## Running instructions

The `phaseinfo` library produces a shared library `phaseinfo.so` that can be dynamically linked to when running ASPECT. To enable linking, just add the following line to the ASPECT input parameter file:

	...
	set Additional shared libraries = /path/to/libphaseinfo.so
	...
	
A basic example implementation of the plugin can be found in `data/aspect-prm-files/example.prm`.

The input parameters enabled by this shared library are explained in more detail [here](https://github.com/cward97/perplex-aspect/wiki/ASPECT-input-parameters).

## Testing

To compile the code with tests enabled, the `PERPLEX_BUILD_TESTING` and `PHASEINFO_BUILD_TESTING` options have to be set to `ON` (default is `OFF`). The tests may then be run with CTest. For example:

	cmake -DAspect_DIR=/path/to/aspect \
	      -DPERPLEX_BUILD_TESTING=ON -DPHASEINFO_BUILD_TESTING=ON \
	      ..
	make -j<N>
	ctest

## Project layout

	data/		ASPECT parameter files and Perple_X data files
	perplex/	Perple_X wrapper
	  extern/	Perple_X source code
	  include/	public header files
	  src/		source code
	  test/		unit tests
	phaseinfo/	ASPECT plugin
	  src/		source code
