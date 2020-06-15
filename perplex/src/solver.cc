//TODO sort these out
#include <algorithm>
#include <cassert>
#include <iterator>

#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <perplex/solver.h>
#include "c_interface.h"
#include "utils.h"


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
}

namespace perplex 
{
  void Solver::init(const std::string filename) 
  {
    // disable stdout to prevent Perple_X dominating stdout
    const int fd = disable_stdout();
    c_interface::solver_init(filename.c_str());
    enable_stdout(fd);

    // initialize private class members
    for (size_t c = 0; c < c_interface::composition_props_get_n_components(); ++c) {
      bulk_composition.push_back(c_interface::bulk_props_get_composition(c));
      composition_names.push_back(std::string{c_interface::composition_props_get_name(c)});
    }

    for (size_t p = 0; p < c_interface::soln_phase_props_get_n(); ++p) {
      phase_names.abbr.push_back(std::string{c_interface::soln_phase_props_get_short_name(p)});
      phase_names.full.push_back(std::string{c_interface::soln_phase_props_get_long_name(p)});
    }

    // TODO
    // Put into a MinimizeResult.init() method.
    minimize_result.phases.resize(c_interface::soln_phase_props_get_n());
    for (size_t p = 0; p < c_interface::soln_phase_props_get_n(); ++p)
      minimize_result.phases[p].composition.resize(c_interface::composition_props_get_n_components());
  }

  const MinimizeResult& Solver::minimize(const double pressure, 
                                         const double temperature)
  {
    // set the initial conditions
    c_interface::solver_set_pressure(utils::convert_pascals_to_bar(pressure));
    c_interface::solver_set_temperature(temperature);
    for (size_t c = 0; c < bulk_composition.size(); ++c)
      c_interface::bulk_props_set_composition(c, bulk_composition[c]);

    // disable stdout to prevent Perple_X dominating stdout
    /* const int fd = disable_stdout(); */
    c_interface::solver_minimize();
    /* enable_stdout(fd); */

    // TODO
    // Put into a MinimizeResult.update() method.
    for (size_t res_phase_idx = 0; res_phase_idx < c_interface::res_phase_props_get_n(); ++res_phase_idx) {
      std::string phase_name{c_interface::res_phase_props_get_name(res_phase_idx)};

        std::vector<std::string>::iterator abbr_it = std::find(phase_names.abbr.begin(), 
	                                                       phase_names.abbr.end(), 
					                       phase_name);
        std::vector<std::string>::iterator full_it = std::find(phase_names.full.begin(), 
	                                phase_names.full.end(), 
					phase_name);

      size_t phase_idx;
      if (abbr_it < phase_names.abbr.end())
	phase_idx = std::distance(phase_names.abbr.begin(), abbr_it);
      else if (full_it < phase_names.full.end())
	phase_idx = std::distance(phase_names.full.begin(), full_it);
      else
	throw std::runtime_error("Phase name not found.");

      Phase &phase { minimize_result.phases[phase_idx] };
      phase.n_moles = c_interface::res_phase_props_get_mol(phase_idx);
      for (size_t c = 0; c < c_interface::composition_props_get_n_components(); ++c) {
	phase.composition[c] = c_interface::res_phase_props_get_composition(phase_idx, c);
      }
    }

    minimize_result.density = c_interface::sys_props_get_density();
    minimize_result.expansivity = c_interface::sys_props_get_expansivity();
    minimize_result.molar_entropy = c_interface::sys_props_get_mol_entropy();
    minimize_result.molar_heat_capacity = c_interface::sys_props_get_mol_heat_capacity();

    return minimize_result;
  }

  const std::vector<double>& Solver::get_bulk_composition() const 
  {
    return bulk_composition;
  }

  void Solver::set_bulk_composition(const std::vector<double> &bulk_composition)
  {
    assert(this->bulk_composition.size() == bulk_composition.size());
    this->bulk_composition = bulk_composition;
  }

  const std::vector<std::string>& Solver::get_composition_names() const {
    return composition_names;
  }

  const PhaseNames& Solver::get_phase_names() const {
    return phase_names;
  }
}
