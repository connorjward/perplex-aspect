#include <iostream>
#include <fstream>
#include <unistd.h>
#include "meemum_wrapper.h"

// constants
constexpr double T0 { 1500 };
constexpr double P0 { 50000 };
constexpr double END_P { 0 };
constexpr double DP { 10000 };

constexpr double g_L { 100 };

// globals
double g_T { T0 }, g_p { P0 };
MeemumWrapper g_meemum;

double bisect(double (*)(double), double, double, double, unsigned);
double compare_melt_frac(double);
double calc_dT(double);
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

	g_T += dT2;

    }
    myfile.close();
}

double calc_dT(double melt_frac) {
    // include latent heat
    // these should be constants...
    const double alpha = g_meemum.expansivity();
    const double Cp = g_meemum.heat_capacity();
    const double rho = g_meemum.density();

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
    g_meemum.minimize(T, p);
    double melt_frac_after { 0.0 };
//   for (int i = 0; i < g_meemum.n_phases(); i++) {
//       if (g_meemum.is_melt(&i))
//           melt_frac_after += g_meemum.phase_amount(&i);
//   }

    return melt_frac_after - melt_change;
}

double compare_melt_frac(const double melt_change_before) {
    return calc_melt_change(melt_change_before) - melt_change_before;
}

