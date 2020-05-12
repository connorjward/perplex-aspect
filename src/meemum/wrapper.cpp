#include <cstddef>
#include <vector>
#include "wrapper.hpp"
#include "calcs.hpp"
#include "props.hpp"

using namespace meemum::props;

MeemumWrapper::MeemumWrapper(const char* filename) { 
    meemum::init(filename);
}

MinimizeResult* MeemumWrapper::minimize(double T, double p) { 
    meemum::minimize(T, p); 

    size_t n { n_phases() };

    Phase** phases { new Phase*[n] };
    
    for (size_t i = 0; i < n; i++) { 
	size_t phase { i + 1 };

	char* name = phase_name(phase);

	phases[i] = new Phase {
	    name,

	    phase_weight_frac(phase),
	    phase_vol_frac(phase),
	    phase_mol_frac(phase),
	    phase_mol(phase),
	};
    }

    return new MinimizeResult {
	sys_density(),
	sys_expansivity(),
	sys_mol_entropy(),
	sys_mol_heat_capacity(),

	phases,
    };
}

std::vector<char*> MeemumWrapper::solution_models() {
    std::vector<char*> models; 

    for (size_t id = 1; id <= n_soln_models(); id++) {
        models.push_back(abbr_soln_name(id));
    }

    return models;
}

