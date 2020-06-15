#include <cstddef>
#include <string>
#include <vector>

namespace perplex 
{

  /**
   * docstring
   */
  struct Phase 
  {
    std::string name;
    double n_moles;
    std::vector<double> composition;
  };

  /**
   * docstring
   */
  struct MinimizeResult 
  {
    double density;
    double expansivity;
    double molar_entropy;
    double molar_heat_capacity;

    std::vector<Phase> phases;
  };

  /**
   * docstring
   */
  class MeemumSolver 
  {
    public:
      /**
       * docstring
       */
      void init(const std::string perplex_filename);

      /**
       * Perform the minimization.
       *
       * @param pressure    The pressure (Pa)
       * @param temperature The temperature (K)
       */
      MinimizeResult minimize(const double pressure, 
	                      const double temperature) const;

      /**
       * docstring
       */
      std::vector<double> get_composition() const;

      /**
       * docstring
       */
      void set_composition(std::vector<double> &composition);

      /**
       * docstring
       */
      std::vector<std::string> get_composition_component_names() const;

      /**
       * docstring
       */
      unsigned int get_n_solution_phases() const;

      /**
       * docstring
       */
      std::vector<std::string> get_solution_phase_names() const;

    private:
      /**
       * docstring
       */
      std::vector<double> composition_;

      /**
       * docstring
       */
      std::vector<std::string> composition_component_names_;

      /**
       * docstring
       */
      std::vector<std::string> solution_phase_names_;
  };

} // end namespace
