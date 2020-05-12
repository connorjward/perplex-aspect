#include <gtest/gtest.h>
#include "meemum/wrapper.hpp"

#include <iostream>

class MeemumWrapperTest : public ::testing::Test {
    protected:
	MinimizeResult *res1, *res2;

	void SetUp() {
	    MeemumWrapper wrapper { "khgp" };
	    res1 = wrapper.minimize(1500, 25000);
	    res2 = wrapper.minimize(1200, 30000);
	}
};

TEST_F(MeemumWrapperTest, CheckPhaseName) {
    EXPECT_STREQ(res1->phases[3]->name, "Opx");
    EXPECT_STREQ(res2->phases[2]->name, "Ol");
}

TEST_F(MeemumWrapperTest, CheckPhaseWeightFrac) {
    EXPECT_NEAR(res1->phases[3]->weight_frac*100, 10.56, 5e-3);
    EXPECT_NEAR(res2->phases[2]->weight_frac*100, 61.45, 5e-3);
}

TEST_F(MeemumWrapperTest, CheckPhaseVolFrac) {
    EXPECT_NEAR(res1->phases[3]->vol_frac*100, 10.68, 5e-3);
    EXPECT_NEAR(res2->phases[2]->vol_frac*100, 62.16, 5e-3);
}

TEST_F(MeemumWrapperTest, CheckPhaseMolFrac) {
    EXPECT_NEAR(res1->phases[3]->mol_frac*100, 8.83, 5e-3);
    EXPECT_NEAR(res2->phases[2]->mol_frac*100, 73.75, 5e-3);
}

/*
TEST(PhasePropsTest, CheckMol) {
    EXPECT_NEAR(res1->phases[3]->mol, 2.62, 5e-3);
    EXPECT_NEAR(res2->phases[2]->mol, 21.5, 0.05);
}

TEST_F(MeemumWrapperTest, CheckDensity) {
    EXPECT_NEAR(res1->sys_density, 3298.4, 0.05);
    EXPECT_NEAR(res2->sys_density, 3356.9, 0.05);
}

TEST_F(MeemumWrapperTest, CheckExpansivity) {
    EXPECT_NEAR(res1->sys_expansivity, 0.37240e-4, 5e-9);
    EXPECT_NEAR(res2->sys_expansivity, 0.34630e-4, 5e-9);
}

TEST(SysPropsTest, CheckMolEntropy) {
    EXPECT_NEAR(res1->sys_mol_entropy, 12536, 0.5);
    EXPECT_NEAR(res2->sys_mol_entropy, 11050, 0.5);
}

TEST(SysPropsTest, CheckMolHeatCapacity) {
    EXPECT_NEAR(res1->sys_mol_heat_capacity, 6496.5, 0.05);
    EXPECT_NEAR(res2->sys_mol_heat_capacity, 6281.6, 0.05);
}
*/
