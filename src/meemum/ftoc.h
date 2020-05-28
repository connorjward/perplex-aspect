#include <cstddef>

extern "C" {
  void c_init(const char*);
  void c_minimize();

  void c_set_pressure(double);
  void c_set_temperature(double);

  double c_get_composition_component(size_t);
  void c_set_composition_component(size_t, double);


  size_t n_soln_models();
  char* abbr_soln_name(size_t);
  char* full_soln_name(size_t);

  size_t n_phases();
  char* phase_name(size_t);
  double phase_weight_frac(size_t);
  double phase_vol_frac(size_t);
  double phase_mol_frac(size_t);
  double phase_mol(size_t);

  double composition_amount(size_t);

  double sys_density();
  double sys_expansivity();
  double sys_mol_entropy();
  double sys_mol_heat_capacity();

}
