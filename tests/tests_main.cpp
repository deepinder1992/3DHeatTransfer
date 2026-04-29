#include <gtest/gtest.h>

int main(int argc, char** argv) {
    std::cout << "=== 3DHeatTransfer Test Suite ===\n\n";
    
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}