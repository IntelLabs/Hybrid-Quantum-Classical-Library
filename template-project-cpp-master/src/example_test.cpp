#include "../third_party/googletest/googletest/include/gtest/gtest.h"

#include "example.hpp"

TEST (AdditionTest, Complete) { 
    EXPECT_EQ(4, addition(2, 2));
    EXPECT_EQ(12, addition(-2, 14));
}