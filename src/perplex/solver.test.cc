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
  EXPECT_NEAR(res.density, 3249.3, 0.05);
  EXPECT_NEAR(res.expansivity, 0.38575e-4, 5e-9);
  EXPECT_NEAR(res.molar_entropy, 11996, 0.5);
  EXPECT_NEAR(res.molar_heat_capacity, 6244.7, 0.05);

  // check solution model info
  EXPECT_EQ(wrapper.solution_phase_names().size(), 4);

  // check phase properties
  EXPECT_STREQ(res.phases[0].name.c_str(), "Cpx(HGP)");
  EXPECT_NEAR(res.phases[0].n_moles, 3.07, 5e-2);
}
