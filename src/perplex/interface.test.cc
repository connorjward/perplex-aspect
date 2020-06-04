#include <gtest/gtest.h>
#include "interface.h"

TEST(FtoCTest, CheckMatchesPerpleXOutput) {
  // initialise
  const char *filename { "test1" };
  const double pressure { 15000 };
  const double temperature { 1000 };
  const double composition[5] { 1.000, 4.000, 25.000, 3.000, 3.500 };

  ftoc::init(filename);
  ftoc::set_pressure(pressure);
  ftoc::set_temperature(temperature);
  for (size_t i = 0; i < 5; i++)
    ftoc::set_composition_component(i, composition[i]);

  // perform minimization
  ftoc::minimize();

  // perform tests
  EXPECT_NEAR(ftoc::get_composition_component(3), 3.000, 5e-4);

  EXPECT_EQ(ftoc::get_n_phases(), 5); 
  EXPECT_STREQ(ftoc::get_phase_name(3), "q");
  EXPECT_NEAR(ftoc::get_phase_weight_frac(3)*100, 43.48, 5e-3);
  EXPECT_NEAR(ftoc::get_phase_mol_frac(3)*100, 70.20, 5e-3);
  EXPECT_NEAR(ftoc::get_phase_mol(3), 17.7, 5e-2);

  EXPECT_NEAR(ftoc::get_sys_density(), 3249.0, 0.05);
  EXPECT_NEAR(ftoc::get_sys_expansivity(), 0.36028e-4, 5e-9);
  EXPECT_NEAR(ftoc::get_sys_mol_entropy(), 4419.1, 0.05);
  EXPECT_NEAR(ftoc::get_sys_mol_heat_capacity(), 2824.5, 0.05);
}

