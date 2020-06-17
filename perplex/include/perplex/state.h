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
       * docstring
       */
      static State& get_instance();

      /**
       *
       */
      void initialize(const std::string& filename);

      /**
       *
       */
      void minimize(const double pressure, const double temperature);

      /**
       * docstring
       */
      size_t get_n_composition_components() const;

      /**
       * docstring
       */
      const std::vector<std::string>& get_composition_component_names() const;

      /**
       * docstring
       */
      size_t get_n_soln_phases() const;

      /**
       * docstring
       */
      const std::vector<std::string>& get_abbr_soln_phase_names() const;

      /**
       * docstring
       */
      const std::vector<std::string>& get_full_soln_phase_names() const;

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
      double get_system_density() const;

      /**
       * docstring
       */
      double get_system_expansivity() const;

      /**
       * docstring
       */
      double get_system_molar_entropy() const;

      /**
       * docstring
       */
      double get_system_molar_heat_capacity() const;

      /**
       * docstring
       */
      size_t get_n_end_phases() const;

      /**
       * docstring
       */
      std::string get_end_phase_name(const size_t end_phase_idx) const;

      /**
       * docstring
       */
      double get_end_phase_weight_frac(const size_t end_phase_idx) const;

      /**
       * docstring
       */
      double get_end_phase_vol_frac(const size_t end_phase_idx) const;

      /**
       * docstring
       */
      double get_end_phase_mol_frac(const size_t end_phase_idx) const;

      /**
       * docstring
       */
      double get_end_phase_mol(const size_t end_phase_idx) const;

      /**
       * docstring
       */
      std::vector<double>
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

      /**
       *
       */
      int get_status() { return status; }

    private:
      /**
       *
       */
      int status;

      /**
       * Construct the class.
       *
       * @remark This constructor is private to enforce the singleton pattern.
       */
      State();
  };
}
