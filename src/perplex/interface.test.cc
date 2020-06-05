#include <gtest/gtest.h>
#include "interface.h"

TEST(InterfaceTest, CheckMatchesPerpleXOutput) {
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
  EXPECT_NEAR(ftoc::get_composition_component(3), 5.880, 5e-4);

  // solution phases
  EXPECT_EQ(ftoc::get_n_soln_models(), 4);
  EXPECT_STREQ(ftoc::get_abbr_soln_name(2), "Ol");
  EXPECT_STREQ(ftoc::get_full_soln_name(3), "Opx(HGP)");

  // check phase information
  EXPECT_EQ(ftoc::get_n_phases(), 3); 
  EXPECT_STREQ(ftoc::get_phase_name(0), "Cpx(HGP)");
  EXPECT_NEAR(ftoc::get_phase_weight_frac(1)*100, 62.02, 5e-3);
  EXPECT_NEAR(ftoc::get_phase_mol_frac(2)*100, 19.70, 5e-3);
  EXPECT_NEAR(ftoc::get_phase_mol(0), 3.07, 5e-2);

  // check system information
  EXPECT_NEAR(ftoc::get_sys_density(), 3249.3, 0.05);
  EXPECT_NEAR(ftoc::get_sys_expansivity(), 0.38575e-4, 5e-9);
  EXPECT_NEAR(ftoc::get_sys_mol_entropy(), 11996, 0.5);
  EXPECT_NEAR(ftoc::get_sys_mol_heat_capacity(), 6244.7, 0.05);
}

