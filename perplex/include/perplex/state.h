#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace perplex
{
  // TODO
  // Make singleton class.
  class State
  {
    public:
      /**
       * docstring
       */
      const size_t n_composition_components;

      /**
       * docstring
       */
      const std::vector<std::string> composition_component_names;

      /**
       * docstring
       */
      const size_t n_soln_phases;

      /**
       * docstring
       */
      const std::vector<std::string> abbr_soln_phase_names;

      /**
       * docstring
       */
      const std::vector<std::string> full_soln_phase_names;

      /**
       * docstring
       */
      State();

      /**
       * docstring
       */
      const std::vector<double>& get_bulk_composition() const;

      /**
       * docstring
       */
      void set_bulk_composition(const std::vector<double>& composition);

      /**
       * docstring
       */
      double get_system_density();

      /**
       * docstring
       */
      double get_system_expansivity();

      /**
       * docstring
       */
      double get_system_molar_entropy();

      /**
       * docstring
       */
      double get_system_molar_heat_capacity();

      /**
       * docstring
       */
      size_t get_n_end_phases();

      /**
       * docstring
       */
      std::string get_end_phase_name(const size_t end_phase_idx);

      /**
       * docstring
       */
      double get_end_phase_weight_frac(const size_t end_phase_idx);

      /**
       * docstring
       */
      double get_end_phase_vol_frac(const size_t end_phase_idx);

      /**
       * docstring
       */
      double get_end_phase_mol_frac(const size_t end_phase_idx);

      /**
       * docstring
       */
      double get_end_phase_mol(const size_t end_phase_idx);

      /**
       * docstring
       */
      const std::vector<double>& 
      get_end_phase_composition(const size_t end_phase_idx) const;

      /**
       * docstring
       */
      const std::string& 
      find_abbr_phase_name(const std::string& end_phase_name) const;

      /**
       * docstring
       */
      const std::string& 
      find_full_phase_name(const std::string& end_phase_name) const;
  };
}
