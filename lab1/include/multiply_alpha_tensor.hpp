#include <vector>

using Matrix = std::vector<std::vector<double>>;

Matrix createMatrix(int rows, int cols, double value = 0);
Matrix multiplyAI(const Matrix& A, const Matrix& B);
