#include <meemum/wrapper.h>
#include "ftoc.h"

void MeemumWrapper::init(const std::string filename) {
  // TODO: set filename before ftoc::init
  ftoc::init(filename.c_str());
}

void MeemumWrapper::minimize(const double pressure, 
                             const double temperature,
			     const std::vector<double> &composition) const {
  // set the temperature and pressure
  ftoc::set_temperature(temperature);
  ftoc::set_pressure(pressure);

  // set the composition
  for (size_t i = 0; i < composition.size(); i++)
    ftoc::set_composition_component(i+1, composition[i]);

  // run the minimization
  ftoc::minimize();
}

size_t MeemumWrapper::n_soln_models() { return ftoc::get_n_soln_models(); }

std::string MeemumWrapper::abbr_soln_name(size_t soln) {
  return std::string { ftoc::get_abbr_soln_name(soln+1) };
}

std::string MeemumWrapper::full_soln_name(size_t soln) {
  return std::string { ftoc::get_full_soln_name(soln+1) };
}

size_t MeemumWrapper::n_phases() { 
  return ftoc::get_n_phases(); 
}

std::string MeemumWrapper::phase_name(const size_t phase) {
  return std::string { ftoc::get_phase_name(phase+1) }; 
}

double MeemumWrapper::phase_weight_frac(const size_t phase) { 
  return ftoc::get_phase_weight_frac(phase+1); 
}

double MeemumWrapper::phase_vol_frac(const size_t phase) { 
  return ftoc::get_phase_vol_frac(phase+1); 
}

double MeemumWrapper::phase_mol_frac(const size_t phase) { 
  return ftoc::get_phase_mol_frac(phase+1); 
}

double MeemumWrapper::phase_mol(const size_t phase) { 
  return ftoc::get_phase_mol(phase+1); 
}

double MeemumWrapper::sys_density() { 
  return ftoc::get_sys_density(); 
}

double MeemumWrapper::sys_expansivity() { return ftoc::get_sys_expansivity(); }

double MeemumWrapper::sys_mol_entropy() { return ftoc::get_sys_mol_entropy(); }

double MeemumWrapper::sys_mol_heat_capacity() { return ftoc::get_sys_mol_heat_capacity(); }

