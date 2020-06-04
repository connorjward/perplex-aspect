#include <gtest/gtest.h>
#include "interface.h"

TEST(FtoCTest, CheckMatchesPerpleXOutput) {
  // initialise
  const char *filename { "test1" };
  const double pressure { 20000 };
  const double temperature { 1500 };
  const double composition[4] { 38.500, 2.820, 50.500, 5.880 };

  ftoc::init(filename);
  ftoc::set_pressure(pressure);
  ftoc::set_temperature(temperature);
  for (size_t i = 0; i < 4; i++)
    ftoc::set_composition_component(i, composition[i]);

  // perform minimization
  ftoc::minimize();

  // perform tests
  EXPECT_NEAR(ftoc::get_composition_component(3), 3.000, 5e-4);

  // solution phases
  EXPECT_EQ(ftoc::get_n_soln_models(), -1);
  EXPECT_STREQ(ftoc::get_abbr_soln_name(2), "???");
  EXPECT_STREQ(ftoc::get_full_soln_name(3), "???");

  // check phase information
  EXPECT_EQ(ftoc::get_n_phases(), 5); 
  EXPECT_STREQ(ftoc::get_phase_name(3), "q");
  EXPECT_NEAR(ftoc::get_phase_weight_frac(3)*100, 43.48, 5e-3);
  EXPECT_NEAR(ftoc::get_phase_mol_frac(3)*100, 70.20, 5e-3);
  EXPECT_NEAR(ftoc::get_phase_mol(3), 17.7, 5e-2);

  // check system information
  EXPECT_NEAR(ftoc::get_sys_density(), 3249.0, 0.05);
  EXPECT_NEAR(ftoc::get_sys_expansivity(), 0.36028e-4, 5e-9);
  EXPECT_NEAR(ftoc::get_sys_mol_entropy(), 4419.1, 0.05);
  EXPECT_NEAR(ftoc::get_sys_mol_heat_capacity(), 2824.5, 0.05);
}

