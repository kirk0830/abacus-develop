#include "../formatter.h"
#include <gtest/gtest.h>

class PhysFmtTest : public testing::Test
{
};

TEST_F(PhysFmtTest, DefaultConstructor) {
    formatter::PhysicalFmt physfmt;
    EXPECT_EQ(physfmt.get_context(), "none");
    EXPECT_NE(physfmt.get_p_formatter(), nullptr);
    EXPECT_FALSE(physfmt.get_decorator_mode());
}

TEST_F(PhysFmtTest, ParameterizedConstructor) {
    formatter::Fmt fmt;
    formatter::PhysicalFmt physfmt("test", &fmt);
    EXPECT_EQ(physfmt.get_context(), "test");
    EXPECT_EQ(physfmt.get_p_formatter(), &fmt);
    EXPECT_TRUE(physfmt.get_decorator_mode());
}

TEST_F(PhysFmtTest, AdjustFormatter) {
    formatter::Fmt fmt;
    formatter::PhysicalFmt physfmt("energy", &fmt);
    physfmt.adjust_formatter();
    EXPECT_EQ(fmt.get_width(), 20);
    EXPECT_EQ(fmt.get_precision(), 10);
    EXPECT_EQ(fmt.get_fillChar(), ' ');
    EXPECT_EQ(fmt.get_fixed(), true);
    EXPECT_EQ(fmt.get_right(), true);
    EXPECT_EQ(fmt.get_error(), false);
    physfmt.adjust_formatter(true);
    EXPECT_EQ(fmt.get_width(), 20);
    EXPECT_EQ(fmt.get_precision(), 10);
    EXPECT_EQ(fmt.get_fillChar(), ' ');
    EXPECT_EQ(fmt.get_fixed(), true);
    EXPECT_EQ(fmt.get_right(), false);
    EXPECT_EQ(fmt.get_error(), false);
}

TEST_F(PhysFmtTest, SetContext) {
    formatter::PhysicalFmt physfmt;
    physfmt.set_context("energy");
    EXPECT_EQ(physfmt.get_context(), "energy");
    EXPECT_EQ(physfmt.get_p_formatter()->get_width(), 20);
    EXPECT_EQ(physfmt.get_p_formatter()->get_precision(), 10);
    EXPECT_EQ(physfmt.get_p_formatter()->get_fillChar(), ' ');
    EXPECT_EQ(physfmt.get_p_formatter()->get_fixed(), true);
    EXPECT_EQ(physfmt.get_p_formatter()->get_right(), true);
    EXPECT_EQ(physfmt.get_p_formatter()->get_error(), false);
}

TEST_F(PhysFmtTest, Setpformatter) {
    formatter::Fmt fmt;
    formatter::PhysicalFmt physfmt;
    physfmt.set_p_formatter(&fmt);
    EXPECT_EQ(physfmt.get_p_formatter(), &fmt);
}

TEST_F(PhysFmtTest, SetDecoratorMode) {
    formatter::PhysicalFmt physfmt;
    physfmt.set_decorator_mode(true);
    EXPECT_TRUE(physfmt.get_decorator_mode());
}
