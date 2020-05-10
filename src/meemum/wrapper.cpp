#include <cstddef>
#include <vector>
#include "meemum/wrapper.hpp"
#include "interface.hpp"

using namespace meemum;

MeemumWrapper::MeemumWrapper(const char* filename) { 
    init(filename);
}

MinimizeResult* MeemumWrapper::minimize(double T, double p) { 
    meemum::minimize(&T, &p); 

    size_t n { n_phases() };

    Phase** phases { new Phase*[n] };
    
    for (size_t i = 0; i < n; i++) { 
	size_t phase { i + 1 };

	char* name { new char[14] };
	load_phase_name(&phase, name);

	phases[i] = new Phase {
	    name,

	    phase_weight_frac(&phase),
	    phase_vol_frac(&phase),
	    phase_mol_frac(&phase),
	    phase_mol(&phase),
	};
    }

    return new MinimizeResult {
	density(),
	entropy(),
	expansivity(),
	heat_capacity(),

	phases,
    };
}

std::vector<char*> MeemumWrapper::solution_models() {
    std::vector<char*> models; 

    for (size_t i = 0; i < n_soln_models(); i++) {
	const size_t id { i + 1 };
	char* name { new char[20] };
	load_abbr_soln_name(&id, name);

	models.push_back(name);
    }

    return models;
}

