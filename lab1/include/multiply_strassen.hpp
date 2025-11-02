#pragma once
#include <cstddef>
#include <vector>

using Matrix = std::vector<std::vector<double>>;

Matrix createMatrix(size_t rows, size_t cols, double value = 0.);
Matrix multiplyStrassen(const Matrix& A, const Matrix& B);
