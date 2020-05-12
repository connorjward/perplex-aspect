#include <gtest/gtest.h>
#include "meemum/wrapper.hpp"

MeemumWrapper* wrapper;

class FooEnvironment : public ::testing::Environment {

    public:
	void SetUp() override {
	    wrapper = new MeemumWrapper("khgp");
	}
};

::testing::Environment* const foo_env =
    ::testing::AddGlobalTestEnvironment(new FooEnvironment);

TEST(MeemumWrapperTest, HGP1) {
    const double T { 1500 };
    const double p { 25000 };

    MinimizeResult* res = wrapper->minimize(T, p);	

    EXPECT_NEAR(res->sys_density, 3298.4, 0.05);
    EXPECT_NEAR(res->sys_mol_entropy, 12536, 0.5);

    Phase* phase = res->phases[3];

    EXPECT_STREQ(phase->name, "Opx");
    EXPECT_NEAR(phase->vol_frac, 0.1068, 5e-5);
    EXPECT_NEAR(phase->mol, 2.62, 5e-3);
}
