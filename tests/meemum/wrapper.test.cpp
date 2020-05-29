#include <gtest/gtest.h>
#include <meemum/wrapper.h>

class MeemumTest : public ::testing::Test {
  protected:
    void SetUp() {
      const std::string filename { "khgp" };
      const double pressure { 25000 };
      const double temperature { 1500 };
      const std::vector<double> composition = { 0.110000e-01, 0.249000, 38.4610, 
	1.77400,
	2.82100,     
	50.5250,     
	5.88200,     
	0.710000E-01 ,
	0.109000     ,
	0.480000E-01 };

      wrapper.init(filename);
      wrapper.minimize(pressure, temperature, composition);

      //res2 = wrapper.minimize(1200, 30000);
    }

    MeemumWrapper wrapper;
};

TEST_F(MeemumTest, CheckPhaseName) {
    EXPECT_STREQ(wrapper.phase_name(3).c_str(), "Opx");
    //EXPECT_STREQ(res2->phases[2]->name, "Ol");
}

TEST_F(MeemumTest, CheckPhaseWeightFrac) {
    EXPECT_NEAR(wrapper.phase_weight_frac(3)*100, 10.56, 5e-3);
    //EXPECT_NEAR(res2->phases[2]->weight_frac*100, 61.45, 5e-3);
}

TEST_F(MeemumTest, CheckPhaseVolFrac) {
    EXPECT_NEAR(wrapper.phase_vol_frac(3)*100, 10.68, 5e-3);
    //EXPECT_NEAR(res2->phases[2]->vol_frac*100, 62.16, 5e-3);
}

TEST_F(MeemumTest, CheckPhaseMolFrac) {
    EXPECT_NEAR(wrapper.phase_mol_frac(3)*100, 8.83, 5e-3);
    //EXPECT_NEAR(res2->phases[2]->mol_frac*100, 73.75, 5e-3);
}

TEST_F(MeemumTest, CheckPhaseMol) {
    EXPECT_NEAR(wrapper.phase_mol(3), 2.62, 5e-3);
    //EXPECT_NEAR(res2->phases[2]->mol, 21.5, 0.05);
}

//TEST_F(MeemumTest, CheckCompositionAmount) {
//    //EXPECT_NEAR(meemum::props::composition_amount(4), 2.62, 5e-3);
//    //EXPECT_NEAR(res2->phases[2]->mol, 21.5, 0.05);
//}

TEST_F(MeemumTest, CheckSysDensity) {
    EXPECT_NEAR(wrapper.sys_density(), 3298.4, 0.05);
    //EXPECT_NEAR(res2->sys_density, 3356.9, 0.05);
}

TEST_F(MeemumTest, CheckSysExpansivity) {
    EXPECT_NEAR(wrapper.sys_expansivity(), 0.37240e-4, 5e-9);
    //EXPECT_NEAR(res2->sys_expansivity, 0.34630e-4, 5e-9);
}

TEST_F(MeemumTest, CheckSysMolEntropy) {
    EXPECT_NEAR(wrapper.sys_mol_entropy(), 12536, 0.5);
    //EXPECT_NEAR(res2->sys_mol_entropy, 11050, 0.5);
}

TEST_F(MeemumTest, CheckSysMolHeatCapacity) {
    EXPECT_NEAR(wrapper.sys_mol_heat_capacity(), 6496.5, 0.05);
    //EXPECT_NEAR(res2->sys_mol_heat_capacity, 6281.6, 0.05);
}
