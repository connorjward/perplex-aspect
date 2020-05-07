#include <iostream>
#include <fstream>
#include <unistd.h>
#include "meemum_wrapper.hpp"

// constants
constexpr double T0 { 1500 };
constexpr double P0 { 50000 };
constexpr double END_P { 0 };
constexpr double DP { 10000 };

constexpr double g_L { 100 };

// globals
double g_T { T0 }, g_p { P0 };

double bisect(double (*)(double), double, double, double, unsigned);
double compare_melt_frac(double);
double calc_dT(double);
double calc_melt_change(double);

int main() {
    std::ofstream myfile;
    myfile.open("mydata.csv");
    myfile << "p,T,S,melt_amount" << std::endl;

    char filenametest[] {"TESTFILENAME"};

    meemum::init(filenametest);
    // needed for extra information eg heat capacity
    meemum::minimize(&g_T, &g_p);


    for ( g_p = P0 ; g_p > END_P; g_p -= DP ) {
	// use bisection to find correct dT
	double min_melt_change = 0.0;
	double max_melt_change = calc_melt_change(min_melt_change);
	const double tol = 0.01;

	//double melt_frac { 0.0 };
	double melt_frac = bisect(&compare_melt_frac, min_melt_change, max_melt_change, tol, 100);

	double dT = calc_dT(melt_frac);

	// update T
	g_T += dT;

        meemum::minimize(&g_T, &g_p);

	// save composition information
	myfile << g_p << "," << g_T << "," 
	    << meemum::get_entropy() << "," << melt_frac << std::endl;
    }
    myfile.close();
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

//  std::cout << alpha << std::endl;
//  std::cout << Cp << std::endl;
//  std::cout << rho << std::endl;
    const double dT = -alpha * g_T * DP *1e5 / Cp / rho - g_L/Cp*melt_frac;
//  std::cout << dT << std::endl;
//  exit(0);

    return dT;
}

double calc_melt_change(double melt_change) {
    double dT = calc_dT(melt_change);
    double T = g_T + dT;
    double p = g_p;
    meemum::minimize(&T, &p);
    double melt_frac_after { 0.0 };
    for (int i = 0; i < meemum::get_n_phases(); i++) {
	if (meemum::is_melt(&i))
	    melt_frac_after += meemum::get_phase_amount(&i);
    }

    return melt_frac_after - melt_change;
}

double compare_melt_frac(const double melt_change_before) {
    return calc_melt_change(melt_change_before) - melt_change_before;
}

