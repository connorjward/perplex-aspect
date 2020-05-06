#include <iostream>
#include <unistd.h>
#include "meemum_wrapper.hpp"

// constants
constexpr double MIN_P {1.0};
constexpr double MAX_P {2.0};
constexpr double DP {0.01};
constexpr double T0 {500.0};
constexpr double g_L {1.0};

// globals
double g_T, g_p;

double bisect(double (*)(double), double, double, double, unsigned);
double compare_melt_frac(double);
double calc_dT(double);
double calc_melt_change(double);

int main() {
    meemum::init();
    g_T = T0;
    // needed for extra information eg heat capacity
    double tmpP = 1.0;
    double tmpT = 500.0;
    meemum::minimize(&tmpT, &tmpP);

    for ( g_p = MIN_P; g_p < MAX_P; g_p += DP ) {
	// save composition information
	std::cout << "some stuff here" << std::endl;

	// use bisection to find correct dT
	double min_melt_change = 0.0;
	double max_melt_change = calc_melt_change(min_melt_change);
	std::cout << max_melt_change << std::endl;
	sleep(1);
	const double tol = 0.01;

	double melt_frac = bisect(&compare_melt_frac, min_melt_change, max_melt_change, tol, 100);

	double dT = calc_dT(melt_frac);

	// update T
	g_T += dT;

	std::cout << g_p << std::endl;
    }
}

double bisect(double (*f)(double), double a, double b, double tol, unsigned int max_iter) {
    double c;

    for ( unsigned int iter = 1; iter < max_iter; iter++ ) {
	// find midpoint
	c = (a + b) / 2;

	if ( (b - a)/2 < tol )
	    return c;

	// new interval
	// if same sign
	if ( f(c) / f(a) > 0 )
	    a = c;
	else
	    b = c;
    }

    // raise exception if unsuccessful
    throw;
}

double calc_dT(double melt_frac) {
    // include latent heat
    // these should be constants...
    const double alpha = meemum::get_expansivity();
    const double Cp = meemum::get_heat_capacity();
    const double rho = meemum::get_density();

    const double dT = -alpha * g_T * DP / Cp / rho - g_L/Cp*melt_frac;
    return dT;
}

double calc_melt_change(double melt_change) {
    double dT = calc_dT(melt_change);
    double T = g_T + dT;
    double p = g_p;
    meemum::minimize(&T, &p);
    const double melt_frac_after = meemum::get_melt_frac();

    return melt_frac_after - melt_change;
}

double compare_melt_frac(const double melt_change_before) {
    return calc_melt_change(melt_change_before) - melt_change_before;
}

