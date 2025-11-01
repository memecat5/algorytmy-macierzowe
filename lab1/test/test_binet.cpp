#include <gtest/gtest.h>
#include "multiply_binet.hpp" 
#include <cmath>
#include <iostream>
#include <random>
#include <iomanip>


// ---------- pomocniczki do testów ----------
Matrix naiveMultiply(const Matrix& A, const Matrix& B) {
    int m = (int)A.size();
    int n = (int)A[0].size();
    int p = (int)B[0].size();
    Matrix C(m, std::vector<double>(p, 0.0));
    for (int i = 0; i < m; ++i)
        for (int k = 0; k < n; ++k)
            for (int j = 0; j < p; ++j)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

bool approxEqualTol(const Matrix& A, const Matrix& B, double eps = 1e-9) {
    if (A.size() != B.size()) return false;
    if (!A.empty() && !B.empty() && A[0].size() != B[0].size()) return false;
    for (size_t i = 0; i < A.size(); ++i)
        for (size_t j = 0; j < A[0].size(); ++j)
            if (std::fabs(A[i][j] - B[i][j]) > eps)
                return false;
    return true;
}

void printMatrixStdout(const Matrix& M, const char* name) {
    std::cout << name << " (" << M.size() << "x" << (M.empty() ? 0 : M[0].size()) << "):\n";
    for (const auto& row : M) {
        for (double v : row) std::cout << std::setw(10) << v << " ";
        std::cout << "\n";
    }
}

// Helper: build random matrix with fixed seed for reproducibility
Matrix randomMatrix(int rows, int cols, std::mt19937& rng,
                    std::uniform_real_distribution<double>& dist) {
    Matrix M(rows, std::vector<double>(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            M[i][j] = dist(rng);
    return M;
}

using Matrix = std::vector<std::vector<double>>;

// Wypisywanie macierzy
void printMatrix(const Matrix& M) {
    for (const auto& row : M) {
        for (int val : row) std::cout << val << " ";
        std::cout << "\n";
    }
}

// ---------------------------------------------------------
// TEST 1 — Nierówne wymiary: 3x2 * 2x4
// ---------------------------------------------------------
TEST(MatrixMultiplicationTest, UnequalSplit_3x2_2x4) {
    Matrix A = {
        {1, 2},
        {3, 4},
        {5, 6}
    };
    Matrix B = {
        {7, 8, 9, 10},
        {11, 12, 13, 14}
    };

    Matrix expected = {
        {29, 32, 35, 38},
        {65, 72, 79, 86},
        {101, 112, 123, 134}
    };

    Matrix result = multiplyMatrix(A, B);
    EXPECT_TRUE(approxEqualTol(result, expected));
}

// ---------------------------------------------------------
// TEST 2 — Nierówne wymiary: 5x3 * 3x2
// ---------------------------------------------------------
TEST(MatrixMultiplicationTest, UnequalSplit_5x3_3x2) {
    Matrix A = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9},
        {10, 11, 12},
        {13, 14, 15}
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
        {22, 23},
        {28, 29}
    };

    Matrix result = multiplyMatrix(A, B);
    EXPECT_TRUE(approxEqualTol(result, expected));
}

// ---------------------------------------------------------
// TEST 3 — Jednowierszowa macierz: 1x5 * 5x1
// ---------------------------------------------------------
TEST(MatrixMultiplicationTest, SingleRowCol_1x5_5x1) {
    Matrix A = { {1, 2, 3, 4, 5} };
    Matrix B = {
        {1},
        {2},
        {3},
        {4},
        {5}
    };

    Matrix expected = { {55} };
    Matrix result = multiplyMatrix(A, B);
    EXPECT_TRUE(approxEqualTol(result, expected));
}

// ---------------------------------------------------------
// TEST 4 — Macierz z liczbami zmiennoprzecinkowymi
// ---------------------------------------------------------
TEST(MatrixMultiplicationTest, FloatingPointValues) {
    Matrix A = {
        {0.5, 1.5},
        {2.5, 3.5}
    };
    Matrix B = {
        {4.0, 5.0},
        {6.0, 7.0}
    };

    Matrix expected = {
        {0.5*4 + 1.5*6, 0.5*5 + 1.5*7},
        {2.5*4 + 3.5*6, 2.5*5 + 3.5*7}
    };

    Matrix result = multiplyMatrix(A, B);
    EXPECT_TRUE(approxEqualTol(result, expected));
}

// ---------------------------------------------------------
// TEST 5 — Nierówne dzielenie: 3x3 * 3x2 (brak potęgi dwójki)
// ---------------------------------------------------------
TEST(MatrixMultiplicationTest, UnequalSplit_3x3_3x2) {
    Matrix A = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };
    Matrix B = {
        {9, 8},
        {7, 6},
        {5, 4}
    };
    Matrix expected = {
        {38, 32},
        {101, 86},
        {164, 140}
    };

    Matrix result = multiplyMatrix(A, B);
    EXPECT_TRUE(approxEqualTol(result, expected));
}

// ---------------------------------------------------------
// TEST 6 — Większa nieregularna macierz (6x5 * 5x3)
// ---------------------------------------------------------
TEST(MatrixMultiplicationTest, LargeUnequal_6x5_5x3) {
    Matrix A = {
        {1, 2, 3, 4, 5},
        {6, 7, 8, 9, 10},
        {11, 12, 13, 14, 15},
        {16, 17, 18, 19, 20},
        {21, 22, 23, 24, 25},
        {26, 27, 28, 29, 30}
    };

    Matrix B = {
        {1, 0, 1},
        {0, 1, 1},
        {1, 1, 0},
        {0, 1, 0},
        {1, 0, 1}
    };

    Matrix expected = {
        {9, 9, 8},
        {24, 24, 23},
        {39, 39, 38},
        {54, 54, 53},
        {69, 69, 68},
        {84, 84, 83}
    };

    Matrix result = multiplyMatrix(A, B);
    EXPECT_TRUE(approxEqualTol(result, expected));
    printMatrix(result);
}

// ---------- TEST: deterministic larger cases (no hardcoded expected) ----------
TEST(RecursiveVsNaive, LargeDeterministic_30x25_25x20) {
    Matrix A = randomMatrix(30, 25, *new std::mt19937(12345), *new std::uniform_real_distribution<double>(-5.0,5.0));
    Matrix B = randomMatrix(25, 20, *new std::mt19937(54321), *new std::uniform_real_distribution<double>(-5.0,5.0));
    Matrix expected = naiveMultiply(A, B);
    Matrix result = multiplyMatrix(A, B);
    if (!approxEqualTol(result, expected, 1e-9)) {
        printMatrixStdout(expected, "Expected");
        printMatrixStdout(result,   "Got     ");
    }
    EXPECT_TRUE(approxEqualTol(result, expected, 1e-9));
}

TEST(RecursiveVsNaive, LargeDeterministic_40x30_30x35) {
    Matrix A = randomMatrix(40, 30, *new std::mt19937(11111), *new std::uniform_real_distribution<double>(-3.0,3.0));
    Matrix B = randomMatrix(30, 35, *new std::mt19937(22222), *new std::uniform_real_distribution<double>(-3.0,3.0));
    Matrix expected = naiveMultiply(A, B);
    Matrix result = multiplyMatrix(A, B);
    if (!approxEqualTol(result, expected, 1e-9)) {
        printMatrixStdout(expected, "Expected");
        printMatrixStdout(result,   "Got     ");
    }
    EXPECT_TRUE(approxEqualTol(result, expected, 1e-9));
}

// ---------- TEST: randomized comparisons (many small/medium cases) ----------
TEST(RecursiveVsNaive, RandomizedComparisons) {
    std::mt19937 rng(20251101); // stały seed => powtarzalne testy
    std::uniform_int_distribution<int> sizeDistSmall(1, 30); // wymiary do 30
    std::uniform_real_distribution<double> valDist(-10.0, 10.0);

    const int TRIES = 12; // liczba losowych przypadków
    for (int t = 0; t < TRIES; ++t) {
        int m = sizeDistSmall(rng);
        int n = sizeDistSmall(rng);
        int p = sizeDistSmall(rng);

        Matrix A = randomMatrix(m, n, rng, valDist);
        Matrix B = randomMatrix(n, p, rng, valDist);

        Matrix expected = naiveMultiply(A, B);
        Matrix result = multiplyMatrix(A, B);

        if (!approxEqualTol(result, expected, 1e-9)) {
            std::cout << "\n---- Random test failed (m=" << m << " n=" << n << " p=" << p << ") ----\n";
            printMatrixStdout(A, "A");
            printMatrixStdout(B, "B");
            printMatrixStdout(expected, "Expected");
            printMatrixStdout(result, "Got");
        }

        EXPECT_TRUE(approxEqualTol(result, expected, 1e-9));
    }
}

// ---------- TEST: a heavier stress test but limited (one case) ----------
TEST(RecursiveVsNaive, StressModerate_60x40_40x30) {
    // UWAGA: test większy — może trwać chwilę. Możesz zmniejszyć rozmiary jeśli za wolno.
    std::mt19937 rng(314159);
    std::uniform_real_distribution<double> valDist(-2.0, 2.0);

    Matrix A = randomMatrix(60, 40, rng, valDist);
    Matrix B = randomMatrix(40, 30, rng, valDist);

    Matrix expected = naiveMultiply(A, B);
    Matrix result = multiplyMatrix(A, B);

    if (!approxEqualTol(result, expected, 1e-9)) {
        std::cout << "\n---- Stress test failed ----\n";
        printMatrixStdout(expected, "Expected");
        printMatrixStdout(result,   "Got     ");
    }

    EXPECT_TRUE(approxEqualTol(result, expected, 1e-9));
}