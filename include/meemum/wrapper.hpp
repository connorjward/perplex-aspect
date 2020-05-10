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
    const double density;
    const double entropy;
    const double expansivity;
    const double heat_capacity;

    Phase** phases;
};

class MeemumWrapper {
    public:
	MeemumWrapper(const char *);
	MinimizeResult* minimize(double, double);
	std::vector<char*> solution_models();
};
