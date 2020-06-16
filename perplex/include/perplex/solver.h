#include <cstddef>
#include <string>
#include <vector>

namespace perplex 
{
  struct CompositionProperties
  {
    std::vector<std::string> component_names;
  }

  struct BulkProperties
  {
    std::vector<double> composition;
  }

  struct SolutionPhaseProperties
  {
    /**
     * docstring
     */
    std::string abbr_name;

    /**
     * docstring
     */
    std::string full_name;
  }

  struct ResultPhaseProperties
  {
    std::string name;
    double n_moles;
    std::vector<double> composition;
  };

  struct SystemPhaseProperties
  {

  }

  /**
   * A data structure containing the result of Solver::minimize().
   */
  struct MinimizeResult
  {
    /**
     * A vector containing the phase information.
     */
    std::vector<PhaseInfo> phases;

    /**
     * The system density (kg/m3).
     */
    double density;

    /**
     * The system expansivity (1/K).
     */
    double expansivity;

    /**
     * The system molar entropy (J/K).
     */
    double molar_entropy;

    /**
     * The system molar heat capacity (J/K).
     */
    double molar_heat_capacity;
  };

  /**
   * A class providing an interface to Perple_X.
   */
  class Solver 
  {
    public:
      /**
       * Initialize the solver.
       *
       * @param filename The Perple_X problem definition file
       */
      void init(const std::string filename);

      /**
       * Perform the minimization.
       *
       * @param pressure    The pressure (Pa).
       * @param temperature The temperature (K).
       * @return            The result of the minimization.
       */
      const MinimizeResult minimize(const double pressure, 
	                             const double temperature);

      /**
       * @return The bulk composition
       */
      const std::vector<double>& get_bulk_composition() const;

      /**
       * @param bulk_composition The new bulk composition
       */
      void set_bulk_composition(const std::vector<double> &bulk_composition);

      /**
       * @return The composition component names.
       */
      const std::vector<std::string>& get_composition_names() const;

      /**
       * @return The phase names.
       */
      const std::vector<PhaseName>& get_solution_phase_names() const;

    private:
      /**
       * A vector storing the bulk composition.
       */
      std::vector<double> bulk_composition;

      std::vector<std::string> composition_names;
      std::vector<PhaseName> solution_phase_names;

      void set_conditions(const double pressure, const double temperature) const;

      void populate_minimize_result(MinimizeResult& minimize_result) const;

      const std::pair<std::string,std::string> 
      get_phase_names(const std::string phase_name) const;
  };
} 
