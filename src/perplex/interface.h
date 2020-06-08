#include <cstddef>

namespace perplex {
  namespace interface {
    extern "C" {
      void init(const char*);
      void minimize();

      void set_pressure(double);
      void set_temperature(double);
      void set_composition_component(size_t, double);

      double get_composition_component(size_t);

      size_t get_n_soln_models();
      char*  get_abbr_soln_name(size_t);
      char*  get_full_soln_name(size_t);

      size_t get_n_phases();
      char*  get_phase_name(size_t);
      double get_phase_weight_frac(size_t);
      double get_phase_vol_frac(size_t);
      double get_phase_mol_frac(size_t);
      double get_phase_mol(size_t);

      double get_sys_density();
      double get_sys_expansivity();
      double get_sys_mol_entropy();
      double get_sys_mol_heat_capacity();
    }
  }
}
