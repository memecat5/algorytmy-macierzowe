#pragma once
#include <vector>

using Matrix = std::vector<std::vector<double>>;

Matrix createMatrix(int rows, int cols, double value = 0.);
Matrix multiplyBinet(const Matrix& A, const Matrix& B);
