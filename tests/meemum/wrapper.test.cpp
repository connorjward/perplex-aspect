#include <gtest/gtest.h>
#include "meemum/calcs.hpp"
#include "meemum/props.hpp"

class MeemumTest : public ::testing::Test {
    protected:
	void SetUp() {
	    meemum::init("khgp");
	    meemum::minimize(1500, 25000);
	    //res2 = wrapper.minimize(1200, 30000);
	}
};

TEST_F(MeemumTest, CheckPhaseName) {
    EXPECT_STREQ(meemum::props::phase_name(4), "Opx");
    //EXPECT_STREQ(res2->phases[2]->name, "Ol");
}

TEST_F(MeemumTest, CheckPhaseWeightFrac) {
    EXPECT_NEAR(meemum::props::phase_weight_frac(4)*100, 10.56, 5e-3);
    //EXPECT_NEAR(res2->phases[2]->weight_frac*100, 61.45, 5e-3);
}

TEST_F(MeemumTest, CheckPhaseVolFrac) {
    EXPECT_NEAR(meemum::props::phase_vol_frac(4)*100, 10.68, 5e-3);
    //EXPECT_NEAR(res2->phases[2]->vol_frac*100, 62.16, 5e-3);
}

TEST_F(MeemumTest, CheckPhaseMolFrac) {
    EXPECT_NEAR(meemum::props::phase_mol_frac(4)*100, 8.83, 5e-3);
    //EXPECT_NEAR(res2->phases[2]->mol_frac*100, 73.75, 5e-3);
}

TEST_F(MeemumTest, CheckPhaseMol) {
    EXPECT_NEAR(meemum::props::phase_mol(4), 2.62, 5e-3);
    //EXPECT_NEAR(res2->phases[2]->mol, 21.5, 0.05);
}

TEST_F(MeemumTest, CheckSysDensity) {
    EXPECT_NEAR(meemum::props::sys_density(), 3298.4, 0.05);
    //EXPECT_NEAR(res2->sys_density, 3356.9, 0.05);
}

TEST_F(MeemumTest, CheckSysExpansivity) {
    EXPECT_NEAR(meemum::props::sys_expansivity(), 0.37240e-4, 5e-9);
    //EXPECT_NEAR(res2->sys_expansivity, 0.34630e-4, 5e-9);
}

TEST_F(MeemumTest, CheckSysMolEntropy) {
    EXPECT_NEAR(meemum::props::sys_mol_entropy(), 12536, 0.5);
    //EXPECT_NEAR(res2->sys_mol_entropy, 11050, 0.5);
}

TEST_F(MeemumTest, CheckSysMolHeatCapacity) {
    EXPECT_NEAR(meemum::props::sys_mol_heat_capacity(), 6496.5, 0.05);
    //EXPECT_NEAR(res2->sys_mol_heat_capacity, 6281.6, 0.05);
}
