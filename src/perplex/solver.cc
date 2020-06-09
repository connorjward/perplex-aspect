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

namespace perplex 
{
  void 
  Solver::init(const std::string filename) 
  {
    // Disable stdout
    const int fd = disable_stdout();
    // Read file
    interface::init(filename.c_str());
    // Re-enable stdout
    enable_stdout(fd);

    // Set composition
    for (unsigned int i = 0; i < interface::get_n_composition_components(); ++i) {
      composition_.push_back(interface::get_composition_component(i));
      composition_component_names_.push_back(std::
	                                    string(interface::
					           get_composition_component_name(i)));
    }

    // Set solution_phase_names
    for (unsigned int i = 0; i < interface::get_n_soln_models(); ++i)
      solution_phase_names_.push_back(std::string(interface::get_full_soln_name(i)));
  }

  MinimizeResult 
  Solver::minimize(const double pressure, 
                   const double temperature) const
  {
    // Set the temperature, pressure and composition
    interface::set_pressure(pressure);
    interface::set_temperature(temperature);
    for (unsigned int i = 0; i < composition_.size(); i++)
      interface::set_composition_component(i, composition_[i]);

    // disable Perple_X output by temporarily disabling stdout
    const int fd = disable_stdout();

    // run the minimization
    interface::minimize();

    // re-enable stdout
    enable_stdout(fd);

    // get phase information
    std::vector<Phase> phases;
    for (unsigned int i = 0; i < interface::get_n_phases(); i++) {
      std::vector<double> phase_composition;
      for (unsigned int j = 0; j < interface::get_n_composition_components(); ++j)
	phase_composition.push_back(interface::get_composition_component(j));

      Phase phase { 
	std::string(interface::get_phase_name(i)),
	interface::get_phase_mol(i),
	phase_composition
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

  std::vector<double> 
  Solver::get_composition() const
  {
    return composition_;   
  }

  void
  Solver::set_composition(std::vector<double> &composition)
  {
    /* AssertThrow(this->composition.size() == composition.size()); */
    this->composition_ = composition;
  }

  std::vector<std::string>
  Solver::get_composition_component_names() const
  {
    return composition_component_names_;
  }

  unsigned int
  Solver::get_n_solution_phases() const
  {
    return solution_phase_names_.size();
  }

  std::vector<std::string> 
  Solver::get_solution_phase_names() const 
  {
    return solution_phase_names_;
  }
}
