#include <cstddef>

extern "C" 
{
  /* ------------------------------------------------------------ */
  /* --------------------- SOLVER FUNCTIONS --------------------- */
  /* ------------------------------------------------------------ */

  /**
   * Initialize the solver.
   *
   * @param filename The name of the Perple_X problem definition file.
   */
  void solver_init(const char* filename);

  /**
   * Perform the minimization. 
   */
  void solver_minimize();

  /**
   * @param pressure The pressure used in the minimization (bar).
   */
  void solver_set_pressure(const double pressure);

  /**
   * @param temperature The temperature used in the minimization (K).
   */
  void solver_set_temperature(const double temperature);

  /* ------------------------------------------------------------ */
  /* ------------------ COMPOSITION PROPERTIES ------------------ */
  /* ------------------------------------------------------------ */

  /**
   * @return Number of composition components.
   */
  size_t composition_props_get_n_components();

  /**
   * @param component_idx Composition component index.
   * @return              Name of a composition component.
   */
  char* composition_props_get_name(size_t component_idx);

  /* ----------------------------------------------------------- */
  /* --------------------- BULK PROPERTIES --------------------- */
  /* ----------------------------------------------------------- */

  /**
   * @param component_idx Composition component index.
   * @return              Amount of component (mol). 
   */
  double bulk_props_get_composition(size_t component_idx);

  /**
   * @param component_idx Composition component index.
   * @param amount        Amount of component (mol).
   */
  void bulk_props_set_composition(size_t component_idx, double amount);

  /* --------------------------------------------------------- */
  /* --------------- SOLUTION PHASE PROPERTIES --------------- */
  /* --------------------------------------------------------- */

  /**
   * @return Number of solution phases.
   */
  size_t soln_phase_props_get_n();

  /**
   * @param soln_phase_idx Solution phase index.
   * @return	       The abbreviated name of the solution phase.
   */
  char* soln_phase_props_get_abbr_name(size_t soln_phase_idx);

  /**
   * @param soln_phase_idx Solution phase index.
   * @return	       The full name of the solution phase.
   */
  char* soln_phase_props_get_full_name(size_t soln_phase_idx);

  /* ----------------------------------------------------------- */
  /* ----------------- RESULT PHASE PROPERTIES ----------------- */
  /* ----------------------------------------------------------- */

  /**
   * @return Number of result phases.
   */
  size_t res_phase_props_get_n();

  /**
   * @param res_phase_idx Result phase index.
   * @return              Result phase name.
   *
   * @remark The returned name can be either the short or long version.
   */
  char* res_phase_props_get_name(size_t res_phase_idx);

  /**
   * @param res_phase_idx Result phase index.
   * @return              Result phase fractional weight.
   */
  double res_phase_props_get_weight_frac(size_t res_phase_idx);

  /**
   * @param res_phase_idx Result phase index.
   * @return              Result phase fractional volume.
   */
  double res_phase_props_get_vol_frac(size_t res_phase_idx);

  /**
   * @param res_phase_idx Result phase index.
   * @return              Result phase fractional molar amount.
   */
  double res_phase_props_get_mol_frac(size_t res_phase_idx);

  /**
   * @param res_phase_idx Result phase index.
   * @return              Result phase molar amount.
   */
  double res_phase_props_get_mol(size_t res_phase_idx);

  /**
   * @param res_phase_idx Result phase index.
   * @param component_idx Composition component index.
   * @return              Amount of component in result phase (mol). 
   */
  double res_phase_props_get_composition(size_t res_phase_idx,
					 size_t component_idx);

  /* ----------------------------------------------------------- */
  /* -------------------- SYSTEM PROPERTIES -------------------- */
  /* ----------------------------------------------------------- */

  /**
   * @return System density (kg/m3).
   */
  double sys_props_get_density();

  /**
   * @return System expansivity (1/K).
   */
  double sys_props_get_expansivity();

  /**
   * @return System molar entropy (J/K).
   */
  double sys_props_get_mol_entropy();

  /**
   * @return System molar heat capacity (J/K).
   */
  double sys_props_get_mol_heat_capacity();
} 
