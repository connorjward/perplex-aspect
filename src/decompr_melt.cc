#include <iostream>
#include <fstream>
#include <unistd.h>
#include "meemum_wrapper.h"

// constants
constexpr double T0 { 1400 };
constexpr double P0 { 50000 };
constexpr double END_P { 30000 };
constexpr double DP { 10000 };

constexpr double g_L { 100 };

// globals
double g_T { T0 }, g_p { P0 };
MeemumWrapper g_meemum;

double bisect(double (*)(double), double, double, double, unsigned);
double compare_melt_frac(double);
double calc_dT();
double calc_melt_change(double);

int main(int argc, char* argv[]) {
    std::ofstream myfile;
    myfile.open("mydata.csv");
    myfile << "p,T,S" << std::endl;

    g_meemum.init(argv[1]);

    for ( g_p = P0 ; g_p > END_P; g_p -= DP ) {
	// save composition information
        g_meemum.minimize(g_T, g_p);
	myfile << g_p << "," << g_T << "," 
	    << g_meemum.entropy() << std::endl;


	double dT = calc_dT();
	g_T +=dT;

    }
    myfile.close();
}

double calc_dT() {
	// determine new temperature step
	double S0 = g_meemum.entropy();

	double dT1 = 0;

	double alpha = g_meemum.expansivity();
	double Cp = g_meemum.heat_capacity();
	double rho = g_meemum.density();

	double dT2 = -alpha * g_T * DP *1e5 / Cp / rho;

	g_meemum.minimize(g_T+dT1, g_p-DP);
	double dS1 = g_meemum.entropy() - S0;

	for (size_t i = 0; i < 3; i++) {
	    g_meemum.minimize(g_T+dT2, g_p-DP);
	    double dS2 = g_meemum.entropy() - S0;

	    double dT3 = (dT1*dS2 - dT2*dS1) / (dS2 - dS1);


	    // update values for next iteration
	    dT1 = dT2;
	    dT2 = dT3;

	    dS1 = dS2;
	}
    return dT2;
}

