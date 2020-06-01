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

void MeemumWrapper::init(const std::string filename) {
  // TODO: set filename before ftoc::init
  std::cout << "starting init" << std::endl;
  // disable stdout
  int bak, new_;
  fflush(stdout);
  bak = dup(1);
  /* new_ = open("/dev/null", O_WRONLY); */
  new_ = open("/dev/null", O_WRONLY);
  dup2(new_, 1);
  // read file
  close(new_);
  ftoc::init(filename.c_str());
  // re-enable stdout
  fflush(stdout);
  dup2(bak, 1);
  close(bak);
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
    ftoc::set_composition_component(i+1, composition[i]);

  std::cout << "starting minimize" << std::endl;
  // disable stdout
  int bak, new_;
  fflush(stdout);
  bak = dup(1);
  new_ = open("/dev/null", O_WRONLY);
  dup2(new_, 1);
  close(new_);

  // run the minimization
  ftoc::minimize();

  // re-enable stdout
  fflush(stdout);
  dup2(bak, 1);
  close(bak);
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

double MeemumWrapper::sys_density() const { 
  return ftoc::get_sys_density(); 
}

double MeemumWrapper::sys_expansivity() const { return ftoc::get_sys_expansivity(); }

double MeemumWrapper::sys_mol_entropy() const { return ftoc::get_sys_mol_entropy(); }

double MeemumWrapper::sys_mol_heat_capacity() const { return ftoc::get_sys_mol_heat_capacity(); }

