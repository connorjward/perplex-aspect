#include <cstddef>
#include <string>
#include <vector>

namespace perplex 
{
  /**
   * A data structure storing both the abbreviated and full names of the
   * possible solution phases.
   */
  struct PhaseNames
  {
    /**
     * docstring
     */
    std::vector<std::string> abbr;

    /**
     * docstring
     */
    std::vector<std::string> full;
  };

  /**
   * A data structure containing phase information.
   */
  struct Phase 
  {
    double n_moles;
    std::vector<double> composition;
  };

  /**
   * A data structure containing the result of Solver::minimize().
   */
  struct MinimizeResult
  {
    /**
     * A vector containing the phase information.
     */
    std::vector<Phase> phases;

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
      const MinimizeResult& minimize(const double pressure, 
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
      const PhaseNames& get_phase_names() const;

    private:
      /**
       * A vector storing the bulk composition.
       */
      std::vector<double> bulk_composition;

      /**
       * A vector storing the names of the different composition components.
       */
      std::vector<std::string> composition_names;

      /**
       * A vector storing the names of the different possible phases in the solution.
       */
      PhaseNames phase_names;

      /**
       * A data structure to store the results of the minimization.
       */
      MinimizeResult minimize_result;
  };
} 
