#include <cstddef>
#include "meemum_wrapper.h"

extern "C" {
    void meemum_init(const char *);
    void meemum_minimize(double *, double *);

    int meemum_n_phases();

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

size_t MeemumWrapper::n_phases() {
    return size_t(meemum_n_phases());
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
