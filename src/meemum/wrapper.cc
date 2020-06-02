#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <meemum/wrapper.h>
#include "ftoc.h"

static const int disable_stdout();
static void enable_stdout(const int);

void MeemumWrapper::init(const std::string filename) {
  // TODO: set filename before ftoc::init
  std::cout << "starting init" << std::endl;
  // disable stdout
  const int stdout_descriptor = disable_stdout();
  // read file
  ftoc::init(filename.c_str());
  // re-enable stdout
  enable_stdout(stdout_descriptor);
  std::cout << "init done" << std::endl;
}

MinimizeResult 
MeemumWrapper::minimize(const double pressure, 
                        const double temperature,
			const std::vector<double> &composition) const {
  // set the temperature and pressure
  ftoc::set_pressure(pressure);
  ftoc::set_temperature(temperature);

  // set the composition
  for (size_t i = 0; i < composition.size(); i++)
    ftoc::set_composition_component(i, composition[i]);

  std::cout << "starting minimize" << std::endl;
  // disable stdout
  const int fd = disable_stdout();

  // run the minimization
  ftoc::minimize();

  // re-enable stdout
  enable_stdout(fd);

  std::cout << "minimize done" << std::endl;

  std::vector<Phase> phases;
  for (size_t i = 0; i < ftoc::get_n_phases(); i++) {
    Phase phase { 
      std::string(ftoc::get_phase_name(i)),
      ftoc::get_phase_mol(i),
    };
    phases.push_back(phase);
  }

  return MinimizeResult {
    ftoc::get_sys_density(),
    ftoc::get_sys_expansivity(),
    ftoc::get_sys_mol_entropy(),
    ftoc::get_sys_mol_heat_capacity(),

    phases,
  };
}

size_t MeemumWrapper::n_soln_models() { return ftoc::get_n_soln_models(); }

std::string MeemumWrapper::abbr_soln_name(size_t soln) {
  return std::string { ftoc::get_abbr_soln_name(soln+1) };
}

std::string MeemumWrapper::full_soln_name(size_t soln) {
  return std::string { ftoc::get_full_soln_name(soln+1) };
}

static const int disable_stdout() {
  // flush stdout
  fflush(stdout);

  // get file descriptors
  const int stdout_descriptor = dup(1);
  const int null_descriptor = open("/dev/null", O_WRONLY);

  // reassign stdout to /dev/null
  dup2(null_descriptor, 1);
  close(null_descriptor);

  return stdout_descriptor;
}

static void enable_stdout(const int stdout_descriptor) {
  // flush stdout
  fflush(stdout);
  
  // reassign descriptor
  dup2(stdout_descriptor, 1);
  close(stdout_descriptor);
}

