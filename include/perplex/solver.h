#include <cstddef>
#include <string>
#include <vector>

namespace perplex 
{
struct Phase 
{
  std::string name;
  double n_moles;
  std::vector<double> composition;
};

struct MinimizeResult 
{
  double density;
  double expansivity;
  double molar_entropy;
  double molar_heat_capacity;

  std::vector<Phase> phases;
};

class Solver 
{
 public:
  void init(const std::string);

  /**
   * Perform the minimization.
   *
   * @param pressure    The pressure (Pa)
   * @param temperature The temperature (K)
   */
  MinimizeResult 
  minimize(const double pressure, 
	   const double temperature) const;
  
  std::vector<double>
  get_composition() const;

  void
  set_composition(std::vector<double> &composition);

  std::vector<std::string>
  get_composition_component_names() const;

  unsigned int
  get_n_solution_phases() const;

  std::vector<std::string> 
  get_solution_phase_names() const;

 private:
  std::vector<double> composition_;

  std::vector<std::string> composition_component_names_;

  std::vector<std::string> solution_phase_names_;
};
}
