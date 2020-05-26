# perplex-aspect

An ASPECT plugin that calculates melt and residue composition using the thermodynamic code Perple_X.
Specifically the code MEEMUM is used that can calculate a plethora of thermodynamic properties at a given P-T-X point.

## Prerequisites

### ASPECT

(work in progress)

### Perple_X

A current version of Perple_X may be obtained from its [website](perplex.ethz.ch/).
Perple_X 6.9.0 is also included in this repo as a zip archive. 

Installation is straightforward to do with the provided makefile. The only requirement of this project is that the PERPLEX_DIR environment variable be set pointing to the installation directory.

### Other programs

- `cmake`
- `gcc`

## Installation instructions

An out of source build is recommended:

	mkdir build
	cd build
	cmake ..
	make

## Testing

Testing is done using [Google Test](https://github.com/google/googletest).
Tests are automatically compiled using `cmake` and can be run using the `ctest` command.
Note that tests are currently very slow as a call to MEEMUM needs to be run for every test.

## Design choices

Previous attempts have been made to combine ASPECT and Perple_X.
However, they had a few problems with their design making the code hard to extend and maintain.

The main difference between this implementation and their is that the initial wrapper code is written in Fortran and not C++.
Although it is possible to integrate the existing Fortran code directly into C++ it requires interfacing directly with obscurely named COMMON blocks (e.g. `cxt22`) and hard-coding some 400-line parameter file.
The code I have written presents a much more C\+\+-friendly interface to the user whilst also being easier to maintain and extend since the parameter file issue may be avoided.

