#include <gtest/gtest.h>
#include "ftoc.h"

TEST(FtoCTest, CheckMatchesPerpleXOutput) {
  // initialise
  const char *filename { "khgp" };
  const double pressure { 25000 };
  const double temperature { 1500 };
  const double composition[10] { 
    0.110000e-01, 
      0.249000, 
      38.4610, 
      1.77400, 
      2.82100, 
      50.5250,     
      5.88200,     
      0.710000e-01,
      0.109000,
      0.480000e-01 
  };

  ftoc::init(filename);
  ftoc::set_pressure(pressure);
  ftoc::set_temperature(temperature);

  for (size_t i = 0; i < 10; i++)
    ftoc::set_composition_component(i, composition[i]);

  ftoc::minimize();

  // perform tests
  EXPECT_NEAR(ftoc::get_composition_component(3), 1.774, 5e-4);

  EXPECT_EQ(ftoc::get_n_phases(), 5); 
  EXPECT_STREQ(ftoc::get_phase_name(3), "Opx");
  EXPECT_NEAR(ftoc::get_phase_weight_frac(3)*100, 10.56, 5e-3);
  EXPECT_NEAR(ftoc::get_phase_mol_frac(3)*100, 8.83, 5e-3);
  EXPECT_NEAR(ftoc::get_phase_mol(3), 2.62, 5e-3);

  EXPECT_NEAR(ftoc::get_sys_density(), 3298.4, 0.05);
  EXPECT_NEAR(ftoc::get_sys_expansivity(), 0.37240e-4, 5e-9);
  EXPECT_NEAR(ftoc::get_sys_mol_entropy(), 12536, 0.5);
  EXPECT_NEAR(ftoc::get_sys_mol_heat_capacity(), 6496.5, 0.05);
}

