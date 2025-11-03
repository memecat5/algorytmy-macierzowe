#include <gtest/gtest.h>
#include <iomanip>
#include <vector>
#include <random>
#include "multiply_alpha_tensor.hpp" 

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

// Klasyczne mnożenie
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

// Porównanie macierzy
bool equalMatrix(const Matrix& A, const Matrix& B, double eps = 1e-6) {
    if (A.size() != B.size() || A[0].size() != B[0].size())
        return false;
    for (size_t i = 0; i < A.size(); ++i)
        for (size_t j = 0; j < A[0].size(); ++j)
            if (std::abs(A[i][j] - B[i][j]) > eps)
                return false;
    return true;
}

void checkEqual(const Matrix& A, const Matrix& B) {
    ASSERT_TRUE(equalMatrix(A, B)) << "Macierze różnią się!";
}

// ---------------------------------------------------------------
// TESTY DLA METODY ALPHATENSOR (4x5 * 5x5)
// ---------------------------------------------------------------

TEST(AlphaTensorTest, Random_4x5x5) {
    Matrix A = randomMatrix(4, 5, 1);
    Matrix B = randomMatrix(5, 5, 2);
    checkEqual(verifyMultiplyClassic(A, B), multiplyAI(A, B));
}

TEST(AlphaTensorTest, IdentityRight) {
    Matrix A = randomMatrix(4, 5, 3);
    Matrix B = {{1,0,0,0,0},
                {0,1,0,0,0},
                {0,0,1,0,0},
                {0,0,0,1,0},
                {0,0,0,0,1}};
    checkEqual(verifyMultiplyClassic(A, B), multiplyAI(A, B));
}

TEST(AlphaTensorTest, IdentityLeft) {
    Matrix A = {{1,0,0,0,0},
                {0,1,0,0,0},
                {0,0,1,0,0},
                {0,0,0,1,0}};
    Matrix B = randomMatrix(5, 5, 4);
    checkEqual(verifyMultiplyClassic(A, B), multiplyAI(A, B));
}

TEST(AlphaTensorTest, AllZeros) {
    Matrix A = std::vector<std::vector<double>>(4, std::vector<double>(5, 0));
    Matrix B = randomMatrix(5, 5, 5);
    checkEqual(verifyMultiplyClassic(A, B), multiplyAI(A, B));
}

TEST(AlphaTensorTest, NegativeValues) {
    Matrix A = randomMatrix(4, 5, 6);
    Matrix B = randomMatrix(5, 5, 7);
    for (auto& row : A)
        for (auto& x : row)
            x *= -1;
    checkEqual(verifyMultiplyClassic(A, B), multiplyAI(A, B));
}

TEST(AlphaTensorTest, StructuredPattern) {
    Matrix A = {
        {1, 2, 3, 4, 5},
        {0, 1, 0, 1, 0},
        {5, 4, 3, 2, 1},
        {-1, 0, 1, 0, -1}
    };
    Matrix B = {
        {2, 0, 1, 0, 3},
        {1, 1, 1, 1, 1},
        {0, 2, 0, 2, 0},
        {3, 0, 3, 0, 3},
        {1, 1, 0, 1, 0}
    };
    checkEqual(verifyMultiplyClassic(A, B), multiplyAI(A, B));
}
