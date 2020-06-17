#include <gtest/gtest.h>

#include <perplex/state.h>
#include <perplex/utils.h>

using namespace perplex;

class StateTest : public ::testing::Test {
  protected:
    void SetUp() override {
      const std::string filename{"test1"};
      const double pressure{utils::convert_bar_to_pascals(20000)};
      const double temperature{1500};
      std::vector<double> composition;
      composition.push_back(38.500);
      composition.push_back(2.820);
      composition.push_back(50.500); 
      composition.push_back(5.880);

      State state{State::get_instance()};

      state.initialize(filename);
      state.set_bulk_composition(composition);
      state.minimize(pressure, temperature);
    }
};

TEST_F(StateTest, CheckNSolutionPhases) 
{
  State state{State::get_instance()};
  EXPECT_EQ(state.get_n_soln_phases(), 4);
}

TEST_F(StateTest, CheckDensity) {
  State state{State::get_instance()};
  EXPECT_NEAR(state.get_system_density(), 3249.3, 0.05);
}

TEST_F(StateTest, CheckExpansivity) {
  State state{State::get_instance()};
  EXPECT_NEAR(state.get_system_expansivity(), 0.38575e-4, 5e-9);
}

TEST_F(StateTest, CheckMolarEntropy) {
  State state{State::get_instance()};
  EXPECT_NEAR(state.get_system_molar_entropy(), 11996, 0.5);
}
 
TEST_F(StateTest, CheckMolarHeatCapacity) {
  State state{State::get_instance()};
  EXPECT_NEAR(state.get_system_molar_heat_capacity(), 6244.7, 0.05);
} 


TEST_F(StateTest, CheckEndPhaseNames) 
{
  State state{State::get_instance()};
  EXPECT_STREQ(state.get_end_phase_name(0).c_str(), "Cpx(HGP)");
}

TEST_F(StateTest, CheckPhaseNMoles) 
{
  State state{State::get_instance()};
  EXPECT_NEAR(state.get_end_phase_mol(0), 3.07, 5e-2);
}
