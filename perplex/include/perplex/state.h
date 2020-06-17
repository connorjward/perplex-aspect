#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace perplex
{
  class State
  {
    public:
      /**
       * Retrieve singleton instance of state. A singleton is used because
       * Perple_X relies heavily on global variables (COMMON blocks) so care
       * must be taken to avoid concurrent access to the resources.
       *
       * @return The singleton instance of the state.
       */
      static State& get_instance();

      /**
       * Initialize Perple_X
       *
       * @param filename The Perple_X problem definition file.
       */
      void initialize(const std::string& filename);

      /**
       * Perform the minimization using MEEMUM.
       *
       * @param pressure    The pressure (Pa).
       * @param temperature The temperature (K).
       */
      void minimize(const double pressure, const double temperature);

      /**
       * @return The number of composition components.
       */
      size_t get_n_composition_components() const;

      /**
       * @return The composition component names.
       */
      const std::vector<std::string>& get_composition_component_names() const;

      /**
       * @return The number of solution phases.
       */
      size_t get_n_soln_phases() const;

      /**
       * @return The abbreviated solution phase names.
       */
      const std::vector<std::string>& get_abbr_soln_phase_names() const;

      /**
       * @return The full solution phase names.
       */
      const std::vector<std::string>& get_full_soln_phase_names() const;

      /**
       * @return The bulk composition.
       */
      const std::vector<double>& get_bulk_composition() const;

      /**
       * @param composition The new bulk composition.
       */
      void set_bulk_composition(const std::vector<double>& composition);

      /**
       * @return The system density (kg/m3).
       */
      double get_system_density() const;

      /**
       * @return The system expansivity (1/K).
       */
      double get_system_expansivity() const;

      /**
       * @return The system molar entropy (J/K).
       */
      double get_system_molar_entropy() const;

      /**
       * @return The system molar heat capacity (J/K).
       */
      double get_system_molar_heat_capacity() const;

      /**
       * @return The number of end phases.
       */
      size_t get_n_end_phases() const;

      /**
       * @param end_phase_idx The index of the end phase.
       * @return	      The end phase name.
       */
      std::string get_end_phase_name(const size_t end_phase_idx) const;

      /**
       * @param end_phase_idx The index of the end phase.
       * @return	      The end phase weight fraction.
       */
      double get_end_phase_weight_frac(const size_t end_phase_idx) const;

      /**
       * @param end_phase_idx The index of the end phase.
       * @return	      The end phase volume fraction.
       */
      double get_end_phase_vol_frac(const size_t end_phase_idx) const;

      /**
       * @param end_phase_idx The index of the end phase.
       * @return	      The end phase molar fraction.
       */
      double get_end_phase_mol_frac(const size_t end_phase_idx) const;

      /**
       * @param end_phase_idx The index of the end phase.
       * @return	      The amount of the end phase (mol).
       */
      double get_end_phase_mol(const size_t end_phase_idx) const;

      /**
       * @param end_phase_idx The index of the end phase.
       * @return	      The composition of the end phase.
       */
      std::vector<double>
      get_end_phase_composition(const size_t end_phase_idx) const;

      /**
       * Return the abbreviated form of the given phase name. This method is required
       * because end phases may be expressed using either the abbreviated or full
       * name of the solution phases.
       *
       * @param phase_name The phase name.
       * @return	   The abbreviated phase name.
       */
      const std::string& 
      find_abbr_phase_name(const std::string& phase_name) const;

      /**
       * Return the full form of the given phase name. This method is required
       * because end phases may be expressed using either the abbreviated or full
       * name of the solution phases.
       *
       * @param phase_name The phase name.
       * @return	   The full phase name.
       */
      const std::string& 
      find_full_phase_name(const std::string& phase_name) const;

    private:
      /**
       * A boolean indicating whether Perple_X has been initialized or not.
       */
      bool initialized = false;

      /**
       * A boolean indicating whether minimize() has been called or not.
       */
      bool minimized = false;

      /**
       * Construct the class.
       *
       * @remark This constructor is private to enforce the singleton pattern.
       */
      State() {};
  };
}
