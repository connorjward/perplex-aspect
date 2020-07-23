# PerpleX-ASPECT

An [ASPECT](https://aspect.geodynamics.org/) plugin that analyses phase compositions using [Perple_X](http://www.perplex.ethz.ch/).
The compositions may be tracked in two ways: using particles and using compositional fields*.

This repository contains the source code for interfacing with ASPECT. The Perple_X interface code can be found [here](https://github.com/cward97/perplex-cpp) (and it is built as part of this project).

\*The compositional field approach is extremely slow and should not be used for actual simulations.

## Prerequisites

- [ASPECT](github.com/geodynamics/aspect)

- CMake version `3.11+`

- A Fortran 2003 compatible Fortran compiler (e.g. `gfortran`)

- A C++ 2017 compatible C++ compiler (e.g. `gcc`)

## Installation instructions

Clone the repository (including submodules):

	$ git clone --recurse-submodules https://github.com/cward97/perplex-aspect.git
	$ cd perplex-aspect

Compile the code with CMake (an out of source build as shown below is recommended):

	$ mkdir build
	$ cd build
	$ cmake -DAspect_DIR=/path/to/aspect/ ..
	$ make -j<N>
	
## Running instructions

The code produces a shared library `libperplexaspect.so` that can be dynamically linked to when running ASPECT. To enable linking, just add the following line to the ASPECT input parameter file:

	...
	set Additional shared libraries = /path/to/libperplexaspect.so
	...
	
The plugin-specific input parameters can then be specified in the rest of the file. 
These are explained in more detail in the source code.

## Cookbooks

Example parameter files may be found in the `cookbooks` directory. The most straightforward way to run them is to execute the following commands:
	
	$ cd build
	$ make setup_cookbooks
	$ cd ../cookbooks
	$ ./aspect parameter-file.prm

## Testing

The tests used in ASPECT are simply used to verify that the output of a simulation remains unchanged.
The tests provided in this repository take this approach and so a test consists of a parameter file and a set of output files whose output may be compared.

Similarly to the cookbooks, the tests can be run as follows:

	$ cd build
	$ make setup_tests
	$ cd ../tests
	$ ./aspect parameter-file.prm

## Project layout

	cookbooks/		example parameter files
	data/perplex		Perple_X data files
	external/perplex-cpp	git submodule containing the Perple_X wrapper code	
	include/		header files
	source/			source code
	tests/			test parameter files
