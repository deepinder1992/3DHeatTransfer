#include <gtest/gtest.h>
#include "../include/linearAlgebra.hpp"

//basic initialization test
TEST(SparseMatrixTest, ConstructorInitializesCorrectly) {
    std::size_t n = 5;
    SparseMatrix mat(n);

    EXPECT_EQ(mat.rows(), n);

    EXPECT_EQ(mat.rowPtr().size(), n + 1);
    for (auto v : mat.rowPtr()) {
        EXPECT_EQ(v, 0);
    }

    // values and colIndex should be empty
    EXPECT_TRUE(mat.values().empty());
    EXPECT_TRUE(mat.colIndex().empty());
}

// Test non-const accessors modify data correctly
TEST(SparseMatrixTest, NonConstAccessorsWork) {
    SparseMatrix mat(3);

    mat.values().push_back(1.5);
    mat.colIndex().push_back(2);
    mat.rowPtr()[1] = 1;

    EXPECT_EQ(mat.values().size(), 1);
    EXPECT_EQ(mat.values()[0], 1.5);

    EXPECT_EQ(mat.colIndex().size(), 1);
    EXPECT_EQ(mat.colIndex()[0], 2);

    EXPECT_EQ(mat.rowPtr()[1], 1);
}

// Test const accessors
TEST(SparseMatrixTest, ConstAccessorsWork) {
    SparseMatrix mat(2);

    mat.values().push_back(3.14);
    mat.colIndex().push_back(1);
    mat.rowPtr()[1] = 1;
    mat.rowPtr()[2] = 1;

    const SparseMatrix& constMat = mat;

    EXPECT_EQ(constMat.rows(), 2);

    EXPECT_EQ(constMat.values()[0], 3.14);
    EXPECT_EQ(constMat.colIndex()[0], 1);
    EXPECT_EQ(constMat.rowPtr()[1], 1);
}

// Test CSR structural consistency (basic sanity)
TEST(SparseMatrixTest, CSRStructureConsistency) {
    SparseMatrix mat(3);

    // Simulate a simple CSR structure:
    // row 0 → 1 element
    // row 1 → 2 elements
    // row 2 → 0 elements
    mat.values() = {10.0, 20.0, 30.0};
    mat.colIndex() = {0, 1, 2};
    mat.rowPtr() = {0, 1, 3, 3};

    EXPECT_EQ(mat.rowPtr().size(), 4);
    EXPECT_EQ(mat.values().size(), 3);
    EXPECT_EQ(mat.colIndex().size(), 3);

    // Check row spans
    EXPECT_EQ(mat.rowPtr()[1] - mat.rowPtr()[0], 1);
    EXPECT_EQ(mat.rowPtr()[2] - mat.rowPtr()[1], 2);
    EXPECT_EQ(mat.rowPtr()[3] - mat.rowPtr()[2], 0);
}