#include <cstddef>
#include "meemum/wrapper.hpp"

extern "C" {
    void meemum_init(const char*);
    void meemum_minimize(const double*, const double*);

    void meemum_load_abbr_soln_name(const size_t*, char*);
    void meemum_load_full_soln_name(const size_t*, char*);
    size_t n_soln_models();

    size_t meemum_n_phases();
    void meemum_load_phase_name(size_t*, char*);

    double meemum_phase_weight_frac(size_t*);
    double meemum_phase_vol_frac(size_t*);
    double meemum_phase_mol_frac(size_t*);
    double meemum_phase_mol(size_t*);

    double meemum_density();
    double meemum_entropy();
    double meemum_expansivity();
    double meemum_heat_capacity();
}

MeemumWrapper::MeemumWrapper(const char* filename) { 
    meemum_init(filename);
}

MinimizeResult* MeemumWrapper::minimize(double T, double p) { 
    meemum_minimize(&T, &p); 

    size_t n_phases { meemum_n_phases() };

    Phase** phases { new Phase*[n_phases] };
    
    for (size_t i = 0; i < n_phases; i++) { 
	size_t phase { i + 1 };

	char* name { new char[14] };
	meemum_load_phase_name(&phase, name);

	phases[i] = new Phase {
	    name,

	    meemum_phase_weight_frac(&phase),
	    meemum_phase_vol_frac(&phase),
	    meemum_phase_mol_frac(&phase),
	    meemum_phase_mol(&phase),
	};
    }

    return new MinimizeResult {
	meemum_density(),
	meemum_entropy(),
	meemum_expansivity(),
	meemum_heat_capacity(),

	phases,
    };
}

char** MeemumWrapper::solution_models() {
    const size_t n_solns = n_soln_models(); 

    char** models { new char*[n_solns] }; 
    for (size_t i = 0; i < n_solns; i++) {
	const size_t id { i + 1 };
	char* name { new char[20] };
	meemum_load_abbr_soln_name(&id, name);
	models[i] = name;
    }

    return models;
}

