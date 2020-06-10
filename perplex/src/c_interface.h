#include <cstddef>

namespace perplex 
{
  namespace c_interface 
  {
    extern "C" 
    {
      void init(const char*);
      void minimize();

      void set_pressure(double);
      void set_temperature(double);
      void set_composition_component(size_t, double);

      size_t get_n_composition_components();
      double get_composition_component(size_t);
      char*  get_composition_component_name(size_t);

      size_t get_n_soln_models();
      char*  get_abbr_soln_name(size_t);
      char*  get_full_soln_name(size_t);

      size_t get_n_phases();
      char*  get_phase_name(size_t);
      double get_phase_weight_frac(size_t);
      double get_phase_vol_frac(size_t);
      double get_phase_mol_frac(size_t);
      double get_phase_mol(size_t);
      double get_phase_composition_component(size_t phase_idx,
	  size_t component_idx);

      double get_sys_density();
      double get_sys_expansivity();
      double get_sys_mol_entropy();
      double get_sys_mol_heat_capacity();
    } 
  } 
} 
