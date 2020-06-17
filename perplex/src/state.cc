#include <cassert>
#include <fcntl.h>
#include <stdexcept>
#include <unistd.h>

#include <perplex/state.h>
#include <perplex/utils.h>
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
}

namespace perplex
{
  State& State::get_instance()
  {
    static State instance;
    return instance;
  }

  void State::initialize(const std::string& filename)
  {
    // disable stdout to prevent Perple_X dominating stdout
    const int fd = disable_stdout();
    solver_init(filename.c_str());
    enable_stdout(fd);

    // save that initialization is complete
    initialized = true;
  }

  void State::minimize(const double pressure, const double temperature)
  {
    solver_set_pressure(utils::convert_pascals_to_bar(pressure));
    solver_set_temperature(temperature);

    // TODO fix this bug
    // disable stdout to prevent Perple_X dominating stdout
    /* const int fd = disable_stdout(); */
    solver_minimize();
    /* enable_stdout(fd); */

    // save that the minimization is complete
    minimized = true;
  }

  size_t State::get_n_composition_components() const
  {
    return composition_props_get_n_components();
  }

  const std::vector<std::string>& 
  State::get_composition_component_names() const
  {
    static std::vector<std::string> names;

    names.resize(get_n_composition_components());
    for (size_t i = 0; i < names.size(); ++i)
      names[i] = std::string{composition_props_get_name(i)};
    return names;
  }

  size_t State::get_n_soln_phases() const
  {
    return soln_phase_props_get_n();
  }

  const std::vector<std::string>& 
  State::get_abbr_soln_phase_names() const
  {
    static std::vector<std::string> names;

    names.resize(get_n_soln_phases());
    for (size_t i = 0; i < names.size(); ++i)
      names[i] = std::string(soln_phase_props_get_abbr_name(i));
    return names;
  }

  const std::vector<std::string>& 
  State::get_full_soln_phase_names() const
  {
    static std::vector<std::string> names;

    names.resize(get_n_soln_phases());
    for (size_t i = 0; i < names.size(); ++i)
      names[i] = std::string(soln_phase_props_get_full_name(i));
    return names;
  }

  const std::vector<double>& State::get_bulk_composition() const
  {
    static std::vector<double> composition;
    
    composition.resize(get_n_composition_components());
    for (size_t i = 0; i < composition.size(); ++i)
      composition[i] = bulk_props_get_composition(i);
    return composition;
  }

  void State::set_bulk_composition(const std::vector<double>& composition)
  {
    if (composition.size() != get_n_composition_components())
      throw std::invalid_argument("Specified bulk composition is the wrong size.");

    for (size_t i = 0; i < composition.size(); ++i)
      bulk_props_set_composition(i, composition[i]);
  }

  double State::get_system_density() const
  {
    return sys_props_get_density();
  }

  double State::get_system_expansivity() const
  {
    return sys_props_get_expansivity();
  }

  double State::get_system_molar_entropy() const
  {
    return sys_props_get_mol_entropy();
  }
  
  double State::get_system_molar_heat_capacity() const
  {
    return sys_props_get_mol_heat_capacity();
  }

  size_t State::get_n_end_phases() const
  {
    return res_phase_props_get_n();
  }

  std::string State::get_end_phase_name(size_t end_phase_idx) const
  {
    return std::string(res_phase_props_get_name(end_phase_idx));
  }

  double State::get_end_phase_weight_frac(size_t end_phase_idx) const
  {
    return res_phase_props_get_weight_frac(end_phase_idx);
  }

  double State::get_end_phase_vol_frac(size_t end_phase_idx) const
  {
    return res_phase_props_get_vol_frac(end_phase_idx);
  }

  double State::get_end_phase_mol_frac(size_t end_phase_idx) const
  {
    return res_phase_props_get_mol_frac(end_phase_idx);
  }

  double State::get_end_phase_mol(size_t end_phase_idx) const
  {
    return res_phase_props_get_mol(end_phase_idx);
  }

  std::vector<double>
  State::get_end_phase_composition(const size_t end_phase_idx) const
  {
    std::vector<double> composition(get_n_composition_components());
    for (size_t i = 0; i < composition.size(); ++i)
      composition[i] = res_phase_props_get_composition(end_phase_idx, i);
    return composition;
  }

  const std::string&
  State::find_abbr_phase_name(const std::string& phase_name) const
  {
    std::vector<std::string> abbr_names{get_abbr_soln_phase_names()};
    std::vector<std::string> full_names{get_full_soln_phase_names()};
    
    assert(abbr_names.size() == get_n_soln_phases());
    assert(full_names.size() == get_n_soln_phases());

    for (size_t i = 0; i < get_n_soln_phases(); ++i) {
      if (phase_name == abbr_names[i] || phase_name == full_names[i])
	return abbr_names[i];
    }
    throw std::invalid_argument("The phase name was not found.");
  }

  const std::string&
  State::find_full_phase_name(const std::string& phase_name) const
  {
    std::vector<std::string> abbr_names{get_abbr_soln_phase_names()};
    std::vector<std::string> full_names{get_full_soln_phase_names()};
    
    assert(abbr_names.size() == get_n_soln_phases());
    assert(full_names.size() == get_n_soln_phases());

    for (size_t i = 0; i < get_n_soln_phases(); ++i) {
      if (phase_name == abbr_names[i] || phase_name == full_names[i])
	return full_names[i];
    }
    throw std::invalid_argument("The phase name was not found.");
  }
}
