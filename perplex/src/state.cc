#include <cassert>

#include <perplex/state.h>
#include "c_interface.h"

namespace
{
  std::vector<std::string> get_composition_component_names()
  {
    std::vector<std::string> names(composition_props_get_n_components());
    
    for (size_t i = 0; i < names.size(); ++i)
      names[i] = std::string(composition_props_get_name(i));
    return names;
  }

  std::vector<std::string> get_abbr_soln_phase_names()
  {
    std::vector<std::string> names(soln_phase_props_get_n());
    for (size_t i = 0; i < names.size(); ++i)
      names[i] = std::string(soln_phase_props_get_abbr_name(i));
    return names;
  }

  std::vector<std::string> get_full_soln_phase_names()
  {
    std::vector<std::string> names(soln_phase_props_get_n());
    for (size_t i = 0; i < names.size(); ++i)
      names[i] = std::string(soln_phase_props_get_full_name(i));
    return names;
  }
}

namespace perplex
{
  State::State() : n_composition_components(composition_props_get_n_components()),
                   composition_component_names(get_composition_component_names()),
		   n_soln_phases(soln_phase_props_get_n()),
		   abbr_soln_phase_names(get_abbr_soln_phase_names()),
		   full_soln_phase_names(get_full_soln_phase_names())
  {}

  const std::vector<double>& State::get_bulk_composition() const
  {
    static std::vector<double> composition(n_composition_components);

    for (size_t i = 0; i < composition.size(); ++i)
      composition[i] = bulk_props_get_composition(i);
    return composition;
  }

  void State::set_bulk_composition(const std::vector<double>& composition)
  {
    assert(composition.size() == n_composition_components);

    for (size_t i = 0; i < composition.size(); ++i)
      bulk_props_set_composition(i, composition[i]);
  }

  double State::get_system_density()
  {
    return sys_props_get_density();
  }

  double State::get_system_expansivity()
  {
    return sys_props_get_expansivity();
  }

  double State::get_system_molar_entropy()
  {
    return sys_props_get_mol_entropy();
  }
  
  double State::get_system_molar_heat_capacity()
  {
    return sys_props_get_mol_heat_capacity();
  }

  size_t State::get_n_end_phases()
  {
    return res_phase_props_get_n();
  }

  std::string State::get_end_phase_name(size_t end_phase_idx)
  {
    return std::string(res_phase_props_get_name(end_phase_idx));
  }

  double State::get_end_phase_weight_frac(size_t end_phase_idx)
  {
    return res_phase_props_get_weight_frac(end_phase_idx);
  }

  double State::get_end_phase_vol_frac(size_t end_phase_idx)
  {
    return res_phase_props_get_vol_frac(end_phase_idx);
  }

  double State::get_end_phase_mol_frac(size_t end_phase_idx)
  {
    return res_phase_props_get_mol_frac(end_phase_idx);
  }

  double State::get_end_phase_mol(size_t end_phase_idx)
  {
    return res_phase_props_get_mol(end_phase_idx);
  }

  const std::vector<double>&
  State::get_end_phase_composition(const size_t end_phase_idx) const
  {
    static std::vector<double> composition(n_composition_components);
    for (size_t i = 0; i < composition.size(); ++i)
      composition[i] = res_phase_props_get_composition(end_phase_idx, i);
    return composition;
  }

  const std::string& 
  State::find_abbr_phase_name(const std::string& end_phase_name) const
  {
    for (size_t i = 0; i < n_soln_phases; ++i) {
      if (end_phase_name == abbr_soln_phase_names[i] || 
	  end_phase_name == full_soln_phase_names[i])
	return abbr_soln_phase_names[i];
    }
    assert(false);
  }

  const std::string&
  State::find_full_phase_name(const std::string& end_phase_name) const
  {
    for (size_t i = 0; i < n_soln_phases; ++i) {
      if (end_phase_name == abbr_soln_phase_names[i] || 
	  end_phase_name == full_soln_phase_names[i])
	return full_soln_phase_names[i];
    }
    assert(false);
  }
}
