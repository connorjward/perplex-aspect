#include <gtest/gtest.h>
#include <meemum/wrapper.h>

TEST(MeemumWrapperTest, CheckMatchesPerpleXOutput) {
  const std::string filename("khgp");
  const double pressure { 25000 };
  const double temperature { 1500 };
  std::vector<double> composition;

  composition.push_back(0.110000e-01);
  composition.push_back(0.249000);
  composition.push_back(38.4610); 
  composition.push_back(1.77400);
  composition.push_back(2.82100);     
  composition.push_back(50.5250);     
  composition.push_back(5.88200);     
  composition.push_back(0.710000e-01);
  composition.push_back(0.109000);
  composition.push_back(0.480000e-01);

  MeemumWrapper wrapper;
  wrapper.init(filename);
  MinimizeResult res = wrapper.minimize(pressure, temperature, composition);

  // check system properties
  EXPECT_NEAR(res.density, 3298.4, 0.05);
  EXPECT_NEAR(res.expansivity, 0.37240e-4, 5e-9);
  EXPECT_NEAR(res.molar_entropy, 12536, 0.5);
  EXPECT_NEAR(res.molar_heat_capacity, 6496.5, 0.05);

  // check phase properties
  Phase phase = res.phases[3];

  EXPECT_STREQ(phase.name.c_str(), "Opx");
  EXPECT_NEAR(phase.n_moles, 2.62, 5e-3);
}
