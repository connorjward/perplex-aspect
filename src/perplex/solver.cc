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

#include <perplex/solver.h>
#include "interface.h"

namespace {
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
}

namespace perplex {
  void Solver::init(const std::string filename) {
    // disable stdout
    const int fd = disable_stdout();
    // read file
    interface::init(filename.c_str());
    // re-enable stdout
    enable_stdout(fd);
  }

  MinimizeResult 
    Solver::minimize(const double pressure, 
	const double temperature,
	const std::vector<double> &composition) const {
      // set the temperature and pressure
      interface::set_pressure(pressure);
      interface::set_temperature(temperature);

      // set the composition
      for (size_t i = 0; i < composition.size(); i++)
	interface::set_composition_component(i, composition[i]);

      // disable Perple_X output by temporarily disabling stdout
      const int fd = disable_stdout();

      // run the minimization
      interface::minimize();

      // re-enable stdout
      enable_stdout(fd);

      // get phase information
      std::vector<Phase> phases;
      for (size_t i = 0; i < interface::get_n_phases(); i++) {
	Phase phase { 
	  std::string(interface::get_phase_name(i)),
	    interface::get_phase_mol(i),
	};
	phases.push_back(phase);
      }

      return MinimizeResult {
	interface::get_sys_density(),
	  interface::get_sys_expansivity(),
	  interface::get_sys_mol_entropy(),
	  interface::get_sys_mol_heat_capacity(),
	  phases
      };
    }

  std::vector<std::string> Solver::solution_phase_names() const {
    std::vector<std::string> names;
    for (size_t i = 0; i < interface::get_n_soln_models(); i++) {
      std::string name { interface::get_full_soln_name(i) };
      names.push_back(name);
    }
    return names;
  }
}
