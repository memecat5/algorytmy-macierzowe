#include <gtest/gtest.h>
#include "multiply_binet.hpp"

TEST(MatrixMultiplicationTest, SmallSquare) {
    Matrix A = {{1, 2}, {3, 4}};
    Matrix B = {{5, 6}, {7, 8}};
    Matrix expected = {{19, 22}, {43, 50}};
    EXPECT_TRUE(equalMatrix(multiplyMatrix(A, B), expected));
}

TEST(MatrixMultiplicationTest, RectangularCase) {
    Matrix A = {{1, 2}, {3, 4}, {5, 6}};
    Matrix B = {{7, 8, 9}, {10, 11, 12}};
    Matrix expected = {
        {27, 30, 33},
        {61, 68, 75},
        {95, 106, 117}
    };
    EXPECT_TRUE(equalMatrix(multiplyMatrix(A, B), expected));
}

TEST(MatrixMultiplicationTest, NonSquareCase) {
    Matrix A = {
        {1, 2, 3, 4, 5},
        {6, 7, 8, 9, 10},
        {11, 12, 13, 14, 15}
    };
    Matrix B = {
        {1, 2},
        {3, 4},
        {5, 6},
        {7, 8},
        {9, 10}
    };
    Matrix expected = {
        {95, 110},
        {220, 260},
        {345, 410}
    };
    EXPECT_TRUE(equalMatrix(multiplyMatrix(A, B), expected));
}

TEST(MatrixMultiplicationTest, SingleElement) {
    Matrix A = {{5}};
    Matrix B = {{7}};
    Matrix expected = {{35}};
    EXPECT_TRUE(equalMatrix(multiplyMatrix(A, B), expected));
}

TEST(MatrixMultiplicationTest, UnevenSplit) {
    Matrix A = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9},
        {10, 11, 12}
    };
    Matrix B = {
        {1, 0},
        {0, 1},
        {1, 1}
    };
    Matrix expected = {
        {4, 5},
        {10, 11},
        {16, 17},
        {22, 23}
    };
    EXPECT_TRUE(equalMatrix(multiplyMatrix(A, B), expected));
}
