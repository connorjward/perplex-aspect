#include <gtest/gtest.h>
#include "interface.h"

using namespace perplex;

class InterfaceTest : public ::testing::Test {
  protected:
    void SetUp() override {
      const char *filename { "test1" };
      const double pressure { 20000 };
      const double temperature { 1500 };
      const double composition[4] { 38.500, 2.820, 50.500, 5.880 };

      interface::init(filename);
      interface::set_pressure(pressure);
      interface::set_temperature(temperature);
      for (size_t i = 0; i < 4; i++)
	interface::set_composition_component(i, composition[i]);

      interface::minimize();
    }
};

TEST_F(InterfaceTest, CheckGetNCompositionComponents) {
  EXPECT_EQ(interface::get_n_composition_components(), 4);
}

TEST_F(InterfaceTest, CheckGetCompositionComponent) {
  EXPECT_NEAR(interface::get_composition_component(3), 5.880, 5e-4);
}

TEST_F(InterfaceTest, CheckGetCompositionComponentName) {
  EXPECT_STREQ(interface::get_composition_component_name(0), "SiO2");
  EXPECT_STREQ(interface::get_composition_component_name(1), "CaO");
  EXPECT_STREQ(interface::get_composition_component_name(2), "MgO");
  EXPECT_STREQ(interface::get_composition_component_name(3), "FeO");
}

TEST_F(InterfaceTest, CheckGetNSolnModels) {
  EXPECT_EQ(interface::get_n_soln_models(), 4);
}

TEST_F(InterfaceTest, CheckGetAbbrSolnName) {
  EXPECT_STREQ(interface::get_abbr_soln_name(2), "Ol");
}

TEST_F(InterfaceTest, CheckGetFullSolnName) {
  EXPECT_STREQ(interface::get_full_soln_name(3), "Opx(HGP)");
}

TEST_F(InterfaceTest, CheckGetNPhases) {
  EXPECT_EQ(interface::get_n_phases(), 3); 
}

TEST_F(InterfaceTest, CheckGetPhaseName) {
  EXPECT_STREQ(interface::get_phase_name(0), "Cpx(HGP)");
}

TEST_F(InterfaceTest, CheckGetPhaseWeightFrac) {
  EXPECT_NEAR(interface::get_phase_weight_frac(1)*100, 62.02, 5e-3);
}

TEST_F(InterfaceTest, CheckGetPhaseMolFrac) {
  EXPECT_NEAR(interface::get_phase_mol_frac(2)*100, 19.70, 5e-3);
}

TEST_F(InterfaceTest, CheckGetPhaseMol) {
  EXPECT_NEAR(interface::get_phase_mol(0), 3.07, 5e-2);
}

TEST_F(InterfaceTest, CheckGetSysDensity) {
  EXPECT_NEAR(interface::get_sys_density(), 3249.3, 0.05);
}

TEST_F(InterfaceTest, CheckGetSysExpansivity) {
  EXPECT_NEAR(interface::get_sys_expansivity(), 0.38575e-4, 5e-9);
}

TEST_F(InterfaceTest, CheckGetSysMolEntropy) {
  EXPECT_NEAR(interface::get_sys_mol_entropy(), 11996, 0.5);
}

TEST_F(InterfaceTest, CheckGetSysMolHeatCapacity) {
  EXPECT_NEAR(interface::get_sys_mol_heat_capacity(), 6244.7, 0.05);
}
