#include <cstddef>

class MeemumWrapper {
    public:
	void init(const char *);
	void minimize(double, double);

	size_t n_phases();
	
	double density();
	double entropy();
	double expansivity();
	double heat_capacity();
};
