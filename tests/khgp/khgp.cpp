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

    EXPECT_NEAR(res->density, 3298.4, 0.05);
    EXPECT_NEAR(res->entropy, 12536, 0.5);
}
