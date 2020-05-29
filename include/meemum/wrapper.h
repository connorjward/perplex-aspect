#include <cstddef>
#include <string>
#include <vector>

class MeemumWrapper {
  public:
    void init(const std::string);

    void minimize(const double pressure, 
	          const double temperature, 
		  const std::vector<double> &composition) const;

    size_t n_soln_models();
    std::string abbr_soln_name(size_t);
    std::string full_soln_name(size_t);

    size_t       n_phases();
    std::string  phase_name(size_t);
    double       phase_weight_frac(size_t);
    double       phase_vol_frac(size_t);
    double       phase_mol_frac(size_t);
    double       phase_mol(size_t);

    double sys_density() const;
    double sys_expansivity() const;
    double sys_mol_entropy() const;
    double sys_mol_heat_capacity() const;
};

