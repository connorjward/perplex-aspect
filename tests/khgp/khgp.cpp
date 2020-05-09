#include <gtest/gtest.h>
#include "meemum_wrapper.h"

class MeemumTest1 : public ::testing::Test {
    char test_file[10] { "test1" };
    double T { 500 };
    double p { 10000 };

    protected:
	void SetUp() override {
	    meemum::init(test_file);
	    meemum::minimize(&T, &p);	
	}
};

TEST_F(MeemumTest1, CorrectNPhases) {
    ASSERT_EQ(meemum::get_n_phases(), 5);
}

TEST_F(MeemumTest1, CorrectDensity) {
    ASSERT_EQ(meemum::get_density(), 4);
}
TEST_F(MeemumTest1, CorrectEntropy) {
    ASSERT_EQ(meemum::get_entropy(), 4);
}
