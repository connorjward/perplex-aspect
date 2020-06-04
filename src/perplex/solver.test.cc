#include <gtest/gtest.h>
#include <perplex/solver.h>

TEST(MeemumWrapperTest, CheckMatchesPerpleXOutput) {
  const std::string filename("test1");
  const double pressure { 20000 };
  const double temperature { 1500 };
  std::vector<double> composition;
  composition.push_back(38.500);
  composition.push_back(2.820);
  composition.push_back(50.500); 
  composition.push_back(5.880);

  MeemumWrapper wrapper;
  wrapper.init(filename);
  MinimizeResult res = wrapper.minimize(pressure, temperature, composition);

  // check system properties
  EXPECT_NEAR(res.density, 3249.0, 0.05);
  EXPECT_NEAR(res.expansivity, 0.36028e-4, 5e-9);
  EXPECT_NEAR(res.molar_entropy, 4419.1, 0.05);
  EXPECT_NEAR(res.molar_heat_capacity, 2824.5, 0.05);

  // check solution model info
  EXPECT_EQ(wrapper.solution_phase_names().size(), -1);

  // check phase properties
  Phase phase = res.phases[3];
  EXPECT_STREQ(phase.name.c_str(), "q");
  EXPECT_NEAR(phase.n_moles, 17.7, 5e-2);
}
