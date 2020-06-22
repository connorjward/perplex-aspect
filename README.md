# PerpleX-ASPECT

This repository contains two libraries: 

- `perplex` : A C++ wrapper for the thermodynamic code [Perple_X](perplex.ethz.ch) (written in Fortran 77).

- `phaseinfo` : A plugin for the geodynamical code [ASPECT](aspect.geodynamics.org) that can analyse phase (mineral) amounts and composition.

## Prerequisites

- [ASPECT](github.com/geodynamics/aspect)

- CMake version `3.11+`

- A Fortran 2003 compatible Fortran compiler (e.g. `gfortran`)

- A C++ 2017 compatible C++ compiler (e.g. `gcc`)

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
	
The plugin-specific input parameters can then be specified in the rest of the file. These are explained in more detail [here](https://github.com/cward97/perplex-aspect/wiki/ASPECT-input-parameters).

## Cookbooks

Example parameter files may be found in the `cookbooks` directory. The most straightforward way to run them is to execute the following commands:
	
	cd build
	make setup_cookbooks
	cd ../cookbooks
	./aspect parameter-file.prm

## Project layout

	cookbooks/	example parameter files
	data/perplex	Perple_X data files
	source/		source code
