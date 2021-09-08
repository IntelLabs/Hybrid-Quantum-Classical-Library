#include "example.hpp"

#include <gtest/gtest.h>

TEST (AdditionTest, Complete) {
    EXPECT_EQ(4, addition(2, 2));
    EXPECT_EQ(12, addition(-2, 14));
}

//int main(int argc, char **argv) {
//    ::testing::InitGoogleTest(&argc, argv);
//    return RUN_ALL_TESTS();
//}
