#include <gtest/gtest.h>
#include <perplex/solver.h>

#include "utils.h"

using namespace perplex;

class SolverTest : public ::testing::Test {
  protected:
    void SetUp() override {
      const std::string filename("test1");
      const double pressure { utils::convert_bar_to_pascals(20000) };
      const double temperature { 1500 };
      std::vector<double> composition;
      composition.push_back(38.500);
      composition.push_back(2.820);
      composition.push_back(50.500); 
      composition.push_back(5.880);

      solver.init(filename);
      solver.set_bulk_composition(composition);
      res = solver.minimize(pressure, temperature);
    }

    Solver solver;
    MinimizeResult res;
};

TEST_F(SolverTest, CheckDensity) {
  EXPECT_NEAR(res.density, 3249.3, 0.05);
}

TEST_F(SolverTest, CheckExpansivity) {
  EXPECT_NEAR(res.expansivity, 0.38575e-4, 5e-9);
}

TEST_F(SolverTest, CheckMolarEntropy) {
  EXPECT_NEAR(res.molar_entropy, 11996, 0.5);
}
 
TEST_F(SolverTest, CheckMolarHeatCapacity) {
  EXPECT_NEAR(res.molar_heat_capacity, 6244.7, 0.05);
} 

TEST_F(SolverTest, CheckNSolutionPhaseNames) {
  EXPECT_EQ(solver.get_phase_names().abbr.size(), 4);
}

/* TEST_F(SolverTest, CheckPhaseNames) { */
/*   EXPECT_STREQ(res.phases[0].name.c_str(), "Cpx"); */
/* } */

TEST_F(SolverTest, CheckPhaseNMoles) {
  EXPECT_NEAR(res.phases[0].n_moles, 3.07, 5e-2);
}
