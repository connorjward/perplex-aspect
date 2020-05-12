#include <cstddef>
#include <vector>

struct Phase {
    const char *name;

    const double weight_frac;
    const double vol_frac;
    const double mol_frac;
    const double mol;
};

struct MinimizeResult {
    const double sys_density;
    const double sys_expansivity;
    const double sys_mol_entropy;
    const double sys_mol_heat_capacity;

    Phase** phases;
};

class MeemumWrapper {
    public:
	MeemumWrapper(const char *);
	MinimizeResult* minimize(double, double);
	std::vector<char*> solution_models();
};
