#include "example/example.h"

#include <cmath>
#include <gtest/gtest.h>

TEST(TestExample, log_binomial_coeff) {
    double d = log_binomial_coeff(10, 3);
    ASSERT_NEAR(log(120.0), d, 1e-6);
}
