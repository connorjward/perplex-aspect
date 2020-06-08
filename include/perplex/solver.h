#include <cstddef>
#include <string>
#include <vector>

namespace perplex {
  struct Phase {
    std::string name;
    double n_moles;
  };

  struct MinimizeResult {
    double density;
    double expansivity;
    double molar_entropy;
    double molar_heat_capacity;

    std::vector<Phase> phases;
  };

  class Solver {
    public:
      void init(const std::string);

      MinimizeResult 
      minimize(const double pressure, 
	       const double temperature, 
	       const std::vector<double> &composition) const;

      std::vector<std::string> solution_phase_names() const;
  };
}
