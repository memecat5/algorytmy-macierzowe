#pragma once
#include <vector>

using Matrix = std::vector<std::vector<double>>;

Matrix createMatrix(int rows, int cols, int value = 0);
Matrix multiplyMatrix(const Matrix& A, const Matrix& B);
bool equalMatrix(const Matrix& A, const Matrix& B);
