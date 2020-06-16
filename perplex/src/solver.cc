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

    for (size_t p = 0; p < c_interface::soln_phase_props_get_n(); ++p)
      solution_phase_names.push_back(PhaseName{
	  std::string{c_interface::soln_phase_props_get_short_name(p)},
	  std::string{c_interface::soln_phase_props_get_long_name(p)}}
      );
  }

  const MinimizeResult Solver::minimize(const double pressure, 
                                        const double temperature)
  {
    set_conditions(pressure, temperature);

    // disable stdout to prevent Perple_X dominating stdout
    /* const int fd = disable_stdout(); */
    c_interface::solver_minimize();
    /* enable_stdout(fd); */

    MinimizeResult res;
    populate_result(res);
    return res;
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

  const std::vector<std::string>& Solver::get_composition_names() const 
  {
    return composition_names;
  }

  const std::vector<std::string>& Solver::get_solution_phase_names() const 
  {
    return solution_phase_names;
  }


  void Solver::set_conditions(const double pressure, const double temperature) const
  {
    c_interface::solver_set_pressure(utils::convert_pascals_to_bar(pressure));
    c_interface::solver_set_temperature(temperature);
    for (size_t c = 0; c < bulk_composition.size(); ++c)
      c_interface::bulk_props_set_composition(c, bulk_composition[c]);
  }

  void Solver::populate_minimize_result(MinimizeResult& minimize_result) const
  {
    for (size_t p = 0; p < c_interface::res_phase_props_get_n(); ++p) {
      PhaseInfo phase_info;

      auto phase_names{get_phase_names(c_interface::res_phase_props_get_name(p))};
      phase_info.abbr_name = phase_names.first;
      phase_info.full_name = phase_names.second();
      phase_info.n_moles = c_interface::res_phase_props_get_mol(p);
      for (size_t c = 0; c < c_interface::composition_props_get_n_components(); ++c) {
	phase_info.composition.push_back(c_interface::res_phase_props_get_composition(p, c));
      }
      minimize_result.phases.push_back(phase_info);
    }

    minimize_result.density = c_interface::sys_props_get_density();
    minimize_result.expansivity = c_interface::sys_props_get_expansivity();
    minimize_result.molar_entropy = c_interface::sys_props_get_mol_entropy();
    minimize_result.molar_heat_capacity = c_interface::sys_props_get_mol_heat_capacity();
  }

  const std::pair<std::string,std::string> 
  Solver::get_phase_names(const std::string phase_name) const 
  {
    for (size_t i = 0; i < c_interface::soln_phase_props_get_n(); ++i) {
      std::string abbr_name{c_interface::soln_phase_props_get_short_name(i)};
      std::string full_name{c_interface::soln_phase_props_get_long_name(i)};

      if (phase_name == abbr_name || phase_name == full_name)
	return std::pair(abbr_name, full_name);
    }
    // TODO throw a suitable error here
  }
}
