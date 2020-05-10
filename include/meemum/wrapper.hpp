#include <cstddef>

const struct Phase {
    char *name;

    double weight_frac;
    double vol_frac;
    double mol_frac;
    double mol;
};

const struct MinimizeResult {
    double density;
    double entropy;
    double expansivity;
    double heat_capacity;

    Phase **phases;
};

class MeemumWrapper {
    public:
	void init(const char *);
	const MinimizeResult* minimize(double, double);
};
