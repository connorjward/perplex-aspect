#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <perplex/meemum_solver.h>
#include "c_interface.h"

namespace 
{

  const int disable_stdout() {
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

  void enable_stdout(const int stdout_descriptor) {
    // flush stdout
    fflush(stdout);

    // reassign descriptor
    dup2(stdout_descriptor, 1);
    close(stdout_descriptor);
  }

  double convert_pascals_to_bar(const double pressure_in_pascals) {
    return pressure_in_pascals / 1e6;
  }

}

namespace perplex 
{
  void 
  MeemumSolver::init(const std::string filename) 
  {
    // Disable stdout
    const int fd = disable_stdout();
    // Read file
    c_interface::init(filename.c_str());
    // Re-enable stdout
    enable_stdout(fd);

    // Set composition
    for (unsigned int i = 0; i < c_interface::get_n_composition_components(); ++i) {
      composition_.push_back(c_interface::get_composition_component(i));
      composition_component_names_.push_back(std::
	  string(c_interface::
	    get_composition_component_name(i)));
    }

    // Set solution_phase_names
    for (unsigned int i = 0; i < c_interface::get_n_soln_models(); ++i)
      solution_phase_names_.push_back(std::string(c_interface::get_full_soln_name(i)));
  }

  MinimizeResult 
  MeemumSolver::minimize(const double pressure, 
	                 const double temperature) const
  {
    // Set the temperature, pressure and composition
    c_interface::set_pressure(convert_pascals_to_bar(pressure));
    c_interface::set_temperature(temperature);
    for (unsigned int i = 0; i < composition_.size(); ++i)
      c_interface::set_composition_component(i, composition_[i]);

    // disable Perple_X output by temporarily disabling stdout
    /* const int fd = disable_stdout(); */

    // run the minimization
    c_interface::minimize();

    // re-enable stdout
    /* enable_stdout(fd); */

    // get phase information
    std::vector<Phase> phases;

    for (unsigned int i = 0; i < solution_phase_names_.size(); ++i) {
      std::string short_name(c_interface::get_abbr_soln_name(i));
      std::string long_name(c_interface::get_full_soln_name(i));
      double molar_amount { 0.0 };
      std::vector<double> phase_composition;

      for (unsigned int j = 0; j < c_interface::get_n_phases(); ++j) {
	// check if solution model present in output
	// since the phase name can sometimes be reported as either the short or long versions
	// both are checked for
	std::cout << short_name << std::endl;
	std::cout << long_name << std::endl;
	std::cout << std::string(c_interface::get_phase_name(j)) << std::endl;

	std::string phase_name{c_interface::get_phase_name(j)};
	if (phase_name == short_name || phase_name == long_name) {
	  molar_amount = c_interface::get_phase_mol(j);

	  for (unsigned int k = 0; k < c_interface::get_n_composition_components(); ++k)
	    phase_composition.push_back(c_interface::get_phase_composition_component(j, k));
	}
      }
      phases.push_back(perplex::Phase{short_name, molar_amount, phase_composition});
    }

    return MinimizeResult {
      c_interface::get_sys_density(),
	c_interface::get_sys_expansivity(),
	c_interface::get_sys_mol_entropy(),
	c_interface::get_sys_mol_heat_capacity(),
	phases
    };
  }

  std::vector<double> 
  MeemumSolver::get_composition() const
  {
    return composition_;   
  }

  void
  MeemumSolver::set_composition(std::vector<double> &composition)
  {
    /* AssertThrow(this->composition.size() == composition.size()); */
    this->composition_ = composition;
  }

  std::vector<std::string>
  MeemumSolver::get_composition_component_names() const
  {
    return composition_component_names_;
  }

  unsigned int
  MeemumSolver::get_n_solution_phases() const
  {
    return solution_phase_names_.size();
  }

  std::vector<std::string> 
  MeemumSolver::get_solution_phase_names() const 
  {
    return solution_phase_names_;
  }

} 
