namespace meemum {
    extern "C" {
	void init(const char*);
	void minimize(const double*, const double*);

	char* abbr_soln_name(const size_t*); // testing
	void load_abbr_soln_name(const size_t*, char*);
	void load_full_soln_name(const size_t*, char*);
	size_t n_soln_models();

	size_t n_phases();
	void load_phase_name(size_t*, char*);

	double phase_weight_frac(size_t*);
	double phase_vol_frac(size_t*);
	double phase_mol_frac(size_t*);
	double phase_mol(size_t*);

	double density();
	double entropy();
	double expansivity();
	double heat_capacity();
    }
}
