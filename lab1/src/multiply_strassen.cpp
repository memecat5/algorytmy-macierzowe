#include <algorithm>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include "multiply_strassen.hpp"


inline size_t rows(const Matrix& A) { return A.size(); }
inline size_t cols(const Matrix& A) { return A.empty() ? 0 : A[0].size(); }

// Tworzy pustą macierz o wymiarach rows x cols
Matrix createMatrix(size_t rows, size_t cols, double value) {
    return Matrix(rows, std::vector<double>(cols, value));
}

/* Dodawanie i odejmowanie dla adaptacyjnego Strassena.
    Przekazujemy referencję do stats, aby zliczać operacje.
*/
Matrix addMatrix(const Matrix& A, const Matrix& B, size_t target_rows, size_t target_cols, StrassenStats& stats){

    Matrix C = createMatrix(target_rows, target_cols);
    
    // Pierwsza pętla: kopiowanie/dodawanie A do zerowej macierzy
    // Operacja += jest liczona jako dodawanie
    for(size_t i = 0; i < std::min(rows(A), target_rows); ++i){
        for(size_t j = 0; j < std::min(cols(A), target_cols); ++j){
            C[i][j] += A[i][j];
            stats.additions++; 
        }
    }
    
    // Druga pętla: właściwe dodawanie B
    for(size_t i = 0; i < std::min(rows(B), target_rows); ++i){
        for(size_t j = 0; j < std::min(cols(B), target_cols); ++j){
            C[i][j] += B[i][j];
            stats.additions++;
        }
    }

    return C;
}

Matrix subtractMatrix(const Matrix& A, const Matrix& B, size_t target_rows, size_t target_cols, StrassenStats& stats){
    Matrix C = createMatrix(target_rows, target_cols);
    
    for(size_t i = 0; i < std::min(rows(A), target_rows); ++i){
        for(size_t j = 0; j < std::min(cols(A), target_cols); ++j){
            C[i][j] += A[i][j];
            stats.additions++;
        }
    }
    
    for(size_t i = 0; i < std::min(rows(B), target_rows); ++i){
        for(size_t j = 0; j < std::min(cols(B), target_cols); ++j){
            C[i][j] -= B[i][j];
            stats.additions++; // Odejmowanie liczymy jako operację addytywną
        }
    }

    return C;
}

// Wycinanie podmacierzy (tylko kopiowanie, bez operacji arytmetycznych)
Matrix slice(const Matrix& A, int start_row, int start_col, int rows, int cols) {
    Matrix sub = createMatrix(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            sub[i][j] = A[start_row + i][start_col + j];
    return sub;
}

// Składanie podmacierzy - tu następuje sumowanie wyników
void paste(Matrix& C, const Matrix& sub, int row, int col, StrassenStats& stats) {
    for (size_t i = 0; i < rows(sub); ++i)
        for (size_t j = 0; j < cols(sub); ++j) {
            C[row + i][col + j] += sub[i][j];
            stats.additions++;
        }
}

// Klasyczne mnożenie dla małych macierzy - zlicza zarówno mnożenia jak i dodawania
Matrix multiplyClassic(const Matrix& A, const Matrix& B, StrassenStats& stats) {
    int m = rows(A), n = cols(A), p = cols(B);
    Matrix C = createMatrix(m, p);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < p; ++j) {
            for (int k = 0; k < n; ++k) {
                C[i][j] += A[i][k] * B[k][j];
                stats.multiplications++; // Mnożenie
                stats.additions++;       // Dodawanie do akumulatora
            }
        }
    }
    return C;
}

Matrix strassenRecursive(const Matrix& A, const Matrix& B, StrassenStats& stats) {
    int m = rows(A);
    int n = cols(A);
    int p = cols(B);

    // Przypadek bazowy
    if (n <= 2) return multiplyClassic(A, B, stats);

    // Dzielimy wymiary na pół
    int m2 = m / 2;
    int n2 = n / 2;
    int p2 = p / 2;

    // Tworzenie podmacierzy (slice nie wykonuje operacji arytmetycznych, więc bez stats)
    Matrix A11 = slice(A, 0, 0, m-m2, n-n2);
    Matrix A12 = slice(A, 0, n-n2, m-m2, n2);
    Matrix A21 = slice(A, m-m2, 0, m2, n-n2);
    Matrix A22 = slice(A, m-m2, n-n2, m2, n2);

    Matrix B11 = slice(B, 0, 0, n-n2, p-p2);
    Matrix B12 = slice(B, 0, p-p2, n-n2, p2);
    Matrix B21 = slice(B, n-n2, 0, n2, p-p2);
    Matrix B22 = slice(B, n-n2, p-p2, n2, p2);

    // 7 iloczynów Strassena - przekazujemy stats do add/subtract i rekurencji
    
    auto A11_A22 = addMatrix(A11, A22, m-m2, n-n2, stats);
    auto B11_B22 = addMatrix(B11, B22, n-n2, p-p2, stats);
    
    Matrix P1 = strassenRecursive(A11_A22, B11_B22, stats);
    
    auto A21_A22 = addMatrix(A21, A22, m2, n-n2, stats);
    
    Matrix P2 = strassenRecursive(A21_A22, B11, stats);
    
    auto B12__B22 = subtractMatrix(B12, B22, n-n2, p2, stats);
    
    Matrix P3 = strassenRecursive(A11, B12__B22, stats);
    
    auto B21__B11 = subtractMatrix(B21, B11, n2, p-p2, stats);
    
    Matrix P4 = strassenRecursive(A22, B21__B11, stats);
    
    auto A11_A12 = addMatrix(A11, A12, m-m2, n2, stats);

    Matrix P5 = strassenRecursive(A11_A12, B22, stats);
    
    auto A21__A11 = subtractMatrix(A21, A11, m-m2, n-n2, stats);
    auto B11_B12 = addMatrix(B11, B12, n-n2, p-p2, stats);
    
    Matrix P6 = strassenRecursive(A21__A11, B11_B12, stats);

    auto A12__A22 = subtractMatrix(A12, A22, m-m2, n2, stats);
    auto B21_B22 = addMatrix(B21, B22, n2, p-p2, stats);

    Matrix P7 = strassenRecursive(A12__A22, B21_B22, stats);

    // Składanie wynikowych bloków - zagnieżdżone wywołania add/sub też muszą mieć stats
    Matrix C11 = 
    addMatrix(
        subtractMatrix(
                addMatrix(P1, P4, m-m2, p-p2, stats),
        P5, m-m2, p-p2, stats),
    P7, m-m2, p-p2, stats);

    Matrix C12 = addMatrix(P3, P5, m-m2, p2, stats);
    Matrix C21 = addMatrix(P2, P4, m2, p-p2, stats);

    Matrix C22 =
    addMatrix(
        subtractMatrix(
            addMatrix(P1, P3, m2, p2, stats),
        P2, m2, p2, stats),
    P6, m2, p2, stats);

    // Wynikowa macierz
    Matrix C = createMatrix(m, p);
    paste(C, C11, 0, 0, stats);
    paste(C, C12, 0, p-p2, stats);
    paste(C, C21, m-m2, 0, stats);
    paste(C, C22, m-m2, p-p2, stats);

    return C;
}

// Główna funkcja mnożąca - wersja ze zbieraniem statystyk
Matrix multiplyStrassen(const Matrix& A, const Matrix& B, StrassenStats& stats) {
    if (B.size() != A[0].size() || A[0].size() != A.size() || B.size() != B[0].size())
        throw std::invalid_argument("Niepoprawne rozmiary macierzy, mnożenie niemożliwe");

    return strassenRecursive(A, B, stats);
}

// Przeciążenie dla zachowania kompatybilności (bez zwracania stats na zewnątrz)
Matrix multiplyStrassen(const Matrix& A, const Matrix& B) {
    StrassenStats stats; // Statystyki są zbierane, ale ignorowane po powrocie
    return multiplyStrassen(A, B, stats);
}