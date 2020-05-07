namespace meemum {
    extern "C" {
	int get_n_components();

	void get_component_amount(int *, double *);
	void get_component_name(int *, char *);

	double get_density();
	double get_entropy();
	double get_expansivity();
	double get_heat_capacity();
	double get_melt_frac();
	int get_n_phases();
	double get_phase_amount(int *);

	bool has_melt();
	bool is_melt(int *);

	void init(char *);
	void minimize(double *, double *);
    }
}

