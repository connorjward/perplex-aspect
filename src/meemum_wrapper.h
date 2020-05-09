#include <cstddef>

struct Phase {
    char *name;

    double weight_frac;
    double vol_frac;
    double mol_frac;
    double mol;
};

class MeemumWrapper {
    public:
	void init(const char *);
	void minimize(double, double);

	struct Phase phase(size_t);
	
	double density();
	double entropy();
	double expansivity();
	double heat_capacity();
};
