#include <cstddef>
#include <string>
#include <vector>

class MeemumWrapper {
  public:
    void init(const std::string);
    void minimize(double, double, const std::vector<double>);

    size_t      n_soln_models();
    std::string abbr_soln_name(size_t);
    std::string full_soln_name(size_t);

    size_t       n_phases();
    std::string  phase_name(size_t);
    double       phase_weight_frac(size_t);
    double       phase_vol_frac(size_t);
    double       phase_mol_frac(size_t);
    double       phase_mol(size_t);

    double sys_density();
    double sys_expansivity();
    double sys_mol_entropy();
    double sys_mol_heat_capacity();
};

