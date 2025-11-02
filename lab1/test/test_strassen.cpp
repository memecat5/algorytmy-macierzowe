#include <gtest/gtest.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <random>
#include "multiply_strassen.hpp"

using Matrix = std::vector<std::vector<double>>;

// Tworzy losową macierz
Matrix randomMatrix(int rows, int cols, int seed = 0) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-5.0, 5.0);
    Matrix M(rows, std::vector<double>(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            M[i][j] = dist(gen);
    return M;
}

// Klasyczne mnożenie (dla porównania)
Matrix verifyMultiplyClassic(const Matrix& A, const Matrix& B) {
    int n = A.size();
    int m = A[0].size();
    int p = B[0].size();
    Matrix C(n, std::vector<double>(p, 0.0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < p; ++j)
            for (int k = 0; k < m; ++k)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

// Pomocnicze wypisywanie macierzy
void printMatrix(const Matrix& M, const std::string& name) {
    std::cout << name << " (" << M.size() << "x" << M[0].size() << "):\n";
    for (auto& row : M) {
        for (double x : row)
            std::cout << std::setw(8) << std::setprecision(4) << x << " ";
        std::cout << "\n";
    }
    std::cout << "\n";
}

// Porównanie macierzy z tolerancją
bool equalMatrix(const Matrix& A, const Matrix& B, double eps = 1e-6) {
    if (A.size() != B.size() || A[0].size() != B[0].size())
        return false;
    for (int i = 0; i < (int)A.size(); ++i)
        for (int j = 0; j < (int)A[0].size(); ++j)
            if (std::abs(A[i][j] - B[i][j]) > eps)
                return false;
    return true;
}

void checkEqual(const Matrix& A, const Matrix& B) {
    if (!equalMatrix(A, B)) {
        printMatrix(A, "Expected");
        printMatrix(B, "Got");
    }
    ASSERT_TRUE(equalMatrix(A, B));
}

//
// -------------------- TESTY --------------------
//

// mała macierz 2x2
TEST(StrassenTest, Small_2x2) {
    Matrix A = {{1, 2}, {3, 4}};
    Matrix B = {{5, 6}, {7, 8}};
    Matrix expected = verifyMultiplyClassic(A, B);
    Matrix result = multiplyStrassen(A, B);
    checkEqual(expected, result);
}

// nieparzysta 3x3
TEST(StrassenTest, Odd_3x3) {
    Matrix A = randomMatrix(3, 3, 1);
    Matrix B = randomMatrix(3, 3, 2);
    checkEqual(verifyMultiplyClassic(A, B), multiplyStrassen(A, B));
}

// większa 5x5
TEST(StrassenTest, Odd_5x5) {
    Matrix A = randomMatrix(5, 5, 3);
    Matrix B = randomMatrix(5, 5, 4);
    checkEqual(verifyMultiplyClassic(A, B), multiplyStrassen(A, B));
}

// większa 6x6
TEST(StrassenTest, Even_6x6) {
    Matrix A = randomMatrix(6, 6, 5);
    Matrix B = randomMatrix(6, 6, 6);
    checkEqual(verifyMultiplyClassic(A, B), multiplyStrassen(A, B));
}

// duża 9x9 (rekurencja głębsza)
TEST(StrassenTest, Large_9x9) {
    Matrix A = randomMatrix(9, 9, 7);
    Matrix B = randomMatrix(9, 9, 8);
    checkEqual(verifyMultiplyClassic(A, B), multiplyStrassen(A, B));
}

// test z wartościami ujemnymi
TEST(StrassenTest, NegativeValues) {
    Matrix A = {{-1, 2, -3}, {4, -5, 6}, {-7, 8, -9}};
    Matrix B = {{9, -8, 7}, {-6, 5, -4}, {3, -2, 1}};
    checkEqual(verifyMultiplyClassic(A, B), multiplyStrassen(A, B));
}

// średnia 10x10
TEST(StrassenTest, Medium_10x10) {
    Matrix A = randomMatrix(10, 10, 9);
    Matrix B = randomMatrix(10, 10, 10);
    checkEqual(verifyMultiplyClassic(A, B), multiplyStrassen(A, B));
}

// większa 16x16 (pełna potęga dwójki)
TEST(StrassenTest, PowerOfTwo_16x16) {
    Matrix A = randomMatrix(16, 16, 11);
    Matrix B = randomMatrix(16, 16, 12);
    checkEqual(verifyMultiplyClassic(A, B), multiplyStrassen(A, B));
}

// nieparzysta 25x25 (rekurencja wielopoziomowa)
TEST(StrassenTest, Odd_25x25) {
    Matrix A = randomMatrix(25, 25, 13);
    Matrix B = randomMatrix(25, 25, 14);
    checkEqual(verifyMultiplyClassic(A, B), multiplyStrassen(A, B));
}

// duża 50x50
TEST(StrassenTest, Large_50x50) {
    Matrix A = randomMatrix(50, 50, 15);
    Matrix B = randomMatrix(50, 50, 16);
    checkEqual(verifyMultiplyClassic(A, B), multiplyStrassen(A, B));
}

// większa nieparzysta 75x75
TEST(StrassenTest, Odd_75x75) {
    Matrix A = randomMatrix(75, 75, 17);
    Matrix B = randomMatrix(75, 75, 18);
    checkEqual(verifyMultiplyClassic(A, B), multiplyStrassen(A, B));
}

// duża 100x100
TEST(StrassenTest, Large_100x100) {
    Matrix A = randomMatrix(100, 100, 19);
    Matrix B = randomMatrix(100, 100, 20);
    checkEqual(verifyMultiplyClassic(A, B), multiplyStrassen(A, B));
}
