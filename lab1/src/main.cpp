#include "multiply_binet.hpp"
#include <iostream>

// Wypisywanie macierzy
void printMatrix(const Matrix& M) {
    for (const auto& row : M) {
        for (int val : row) std::cout << val << " ";
        std::cout << "\n";
    }
}

// --- main ---
int main() {
    Matrix A = {
        {2,1,3},
        {7,2,1},
        {3,7,2}
    };
    Matrix B = {
        {1},
        {2},
        {3}
    };


    Matrix C = multiplyBinet(A, B);

    std::cout << "\nWynik A * B:\n";
    printMatrix(C);
}
