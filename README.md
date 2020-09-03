# PerpleX-ASPECT

A set of [ASPECT](https://aspect.geodynamics.org/) plugins that analyse phase composition information using [Perple_X](http://www.perplex.ethz.ch/).
The compositions may be tracked in two ways: using particles and using compositional fields.

This repository contains the source code for interfacing with ASPECT and is built as a shared library.
The Perple_X interface code can be found [here](https://github.com/cward97/perplex-cpp) (and it is built as part of this project).

## Prerequisites

- [ASPECT](github.com/geodynamics/aspect)

- CMake version `2.8.12`

- A Fortran 2003 compatible Fortran compiler (e.g. `gfortran`)

- A C++ 2011 compatible C++ compiler (e.g. `gcc`)


## Installation instructions

Clone the repository (including submodules):

	$ git clone --recurse-submodules https://github.com/cward97/perplex-aspect.git
	$ cd perplex-aspect

Compile the code with CMake (an out of source build as shown below is recommended):

	$ mkdir build
	$ cd build
	$ cmake -DAspect_DIR=/path/to/aspect/ ..
	$ make -j<N>

Note that the same compilers that installed ASPECT and its dependencies must be used here.
	

## Running instructions

The code produces a shared library `libperplexaspect.so` that can be dynamically linked to when running ASPECT. 
To enable linking, just add the following line to the ASPECT input parameter file:

	...
	set Additional shared libraries = /path/to/libperplexaspect.so
	...
	
The plugin-specific input parameters can then be specified in the rest of the file. 
These are explained in more detail in the source code.


## Cookbooks

Example parameter files may be found in the `cookbooks` directory. 
The most straightforward way to run them is to execute the following commands:
	
	$ cd build
	$ make setup_cookbooks
	$ cd ../cookbooks
	$ ./aspect parameter-file.prm


## Testing

The tests used in ASPECT are simply used to verify that the output of a simulation does not change.
Each test used here consists of a parameter file and a set of output files whose output may be compared to future runs of the simulation.

Similarly to the cookbooks, the tests can be run as follows:

	$ cd build
	$ make setup_tests
	$ cd ../tests
	$ ./aspect parameter-file.prm


## Perple_X data sets

Two Perple_X data sets are provided in the repository, both modelling KLB-1 peridotite.
The `klb-1` data set was taken from [here](http://www.perplex.ethz.ch/perplex/examples/example_holland_et_al_2018_melt_model/).
The `simple` data set is practically the same but has been altered to reduce the complexity and run time.


## Useful links

- [Report](https://github.com/cward97/miscada-report) discussing the code (2020).


## Project layout

	cookbooks/		example parameter files
	data/aspect		ASPECT data files
	data/perplex		Perple_X data files
	external/perplex-cpp	Perple_X wrapper (git submodule)	
	include/		header files
	source/			source code
	tests/			test parameter files
