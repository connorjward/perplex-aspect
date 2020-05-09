#include <cstddef>
#include "meemum_wrapper.h"

extern "C" {
    void meemum_init(const char*);
    void meemum_minimize(const double*, const double*);

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

void MeemumWrapper::init(const char* filename) { 
    meemum_init(filename);
}

void MeemumWrapper::minimize(double T, double p) { 
    meemum_minimize(&T, &p); 
}

struct Phase MeemumWrapper::phase(size_t phase_id) {
    char* name { new char[14] };
    meemum_load_phase_name(&phase_id, name);

    const double wt = meemum_phase_weight_frac(&phase_id);
    const double vol = meemum_phase_vol_frac(&phase_id);
    const double mol_frac = meemum_phase_mol_frac(&phase_id);
    const double mol = meemum_phase_mol(&phase_id);

    struct Phase phase { name, wt, vol, mol_frac, mol };

    return phase;
}

double MeemumWrapper::density() {
    return meemum_density();
}

double MeemumWrapper::entropy() {
    return meemum_entropy();
}
double MeemumWrapper::expansivity() {
    return meemum_expansivity();
}

double MeemumWrapper::heat_capacity() {
    return meemum_heat_capacity();
}
