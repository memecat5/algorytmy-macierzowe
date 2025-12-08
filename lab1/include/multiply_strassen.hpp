#pragma once
#include <cstddef>
#include <vector>

using Matrix = std::vector<std::vector<double>>;

// Struktura do przechowywania statystyk
struct StrassenStats {
    long long additions = 0;       // Licznik dodawań (i odejmowań)
    long long multiplications = 0; // Licznik mnożeń
};

Matrix createMatrix(size_t rows, size_t cols, double value = 0.);
Matrix multiplyStrassen(const Matrix& A, const Matrix& B);
Matrix multiplyStrassen(const Matrix& A, const Matrix& B, StrassenStats& stats);
