#include "../../src/domain/rect_domain.h"
#include "../../src/initial_condition/initial_condition.h"
#include "../../src/utils.h"
#include <iostream>
#include <complex>
#include <math.h>
#include "gtest/gtest.h"

TEST(InitialConditionTest, RectangularTest)
{
    auto domain = RectangularDomain(21, 21, 0, 10, 11, -10, 10, -10, 10);
    auto initial_cond_function = [](float x, float y)
    { return std::complex<float>{float(1.)}; }; //  expf(-(x * x + y * y) / (9))

    auto *initial_condition = new InitialCondition(initial_cond_function);
    initial_condition->assign_to_domain(&domain);
    ASSERT_FLOAT_EQ(domain.at(10, 10, 0)->x, 0.);
    ASSERT_FLOAT_EQ(domain.at(10, 10, 0)->y, 0.);
    ASSERT_FLOAT_EQ(domain.at(10, 10, 0)->value.real(), 0.047619049);
    ASSERT_FLOAT_EQ(domain.at(10, 10, 0)->value.imag(), 0.);
}
