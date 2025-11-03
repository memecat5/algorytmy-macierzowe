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

/* 
    Dodawanie i odejmowanie dla adaptacyjnego Strassena, które rozszerzają/ucinają wynik,
    w zależności od tego jaki ma być rozmiar wyniku (target_rows, target_cols)
*/
Matrix addMatrix(const Matrix& A, const Matrix& B, size_t target_rows, size_t target_cols){

    Matrix C = createMatrix(target_rows, target_cols);
    for(size_t i = 0; i < std::min(rows(A), target_rows); ++i){
        for(size_t j = 0; j < std::min(cols(A), target_cols); ++j){
            C[i][j] += A[i][j];
        }
    }
    for(size_t i = 0; i < std::min(rows(B), target_rows); ++i){
        for(size_t j = 0; j < std::min(cols(B), target_cols); ++j){
            C[i][j] += B[i][j];
        }
    }

    return C;
}
Matrix subtractMatrix(const Matrix& A, const Matrix& B, size_t target_rows, size_t target_cols){
        Matrix C = createMatrix(target_rows, target_cols);
    for(size_t i = 0; i < std::min(rows(A), target_rows); ++i){
        for(size_t j = 0; j < std::min(cols(A), target_cols); ++j){
            C[i][j] += A[i][j];
        }
    }
    for(size_t i = 0; i < std::min(rows(B), target_rows); ++i){
        for(size_t j = 0; j < std::min(cols(B), target_cols); ++j){
            C[i][j] -= B[i][j];
        }
    }

    return C;
}

// Wycinanie podmacierzy
Matrix slice(const Matrix& A, int start_row, int start_col, int rows, int cols) {
    Matrix sub = createMatrix(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            sub[i][j] = A[start_row + i][start_col + j];
    return sub;
}

// Składanie podmacierzy
void paste(Matrix& C, const Matrix& sub, int row, int col) {
    for (size_t i = 0; i < rows(sub); ++i)
        for (size_t j = 0; j < cols(sub); ++j)
            C[row + i][col + j] += sub[i][j];
}

// Klasyczne mnożenie dla małych macierzy
Matrix multiplyClassic(const Matrix& A, const Matrix& B) {
    int m = rows(A), n = cols(A), p = cols(B);
    Matrix C = createMatrix(m, p);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < p; ++j)
            for (int k = 0; k < n; ++k)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

Matrix strassenRecursive(const Matrix& A, const Matrix& B) {
    int m = rows(A);
    int n = cols(A);
    int p = cols(B);

    // przypadek bazowy
    if (n <= 2) return multiplyClassic(A, B);


    // Dzielimy każdy wymiar na pół (dzielenie intów, więc tu może wyjść 0)
    int m2 = m / 2;
    int n2 = n / 2;
    int p2 = p / 2;

    /*  Tym razem inny podział         
              n-n2  n2
               ^     ^
        m-m2{ A11 | A12
    A =       ---------
        m2{   A21 | A22
    
              p-p2  p2
               ^     ^
        n-n2{ B11 | B12
    B =       ---------
        n2{   B21 | B22
    */

    Matrix A11 = slice(A, 0, 0, m-m2, n-n2);
    Matrix A12 = slice(A, 0, n-n2, m-m2, n2);
    Matrix A21 = slice(A, m-m2, 0, m2, n-n2);
    Matrix A22 = slice(A, m-m2, n-n2, m2, n2);

    Matrix B11 = slice(B, 0, 0, n-n2, p-p2);
    Matrix B12 = slice(B, 0, p-p2, n-n2, p2);
    Matrix B21 = slice(B, n-n2, 0, n2, p-p2);
    Matrix B22 = slice(B, n-n2, p-p2, n2, p2);

    // 7 iloczynów Strassena

    auto A11_A22 = addMatrix(A11, A22, m-m2, n-n2);
    auto B11_B22 = addMatrix(B11, B22, n-n2, p-p2);
    
    Matrix P1 = strassenRecursive(A11_A22, B11_B22);
    
    auto A21_A22 = addMatrix(A21, A22, m2, n-n2);
    
    Matrix P2 = strassenRecursive(A21_A22, B11);
    
    auto B12__B22 = subtractMatrix(B12, B22, n-n2, p2);
    
    Matrix P3 = strassenRecursive(A11, B12__B22);
    
    auto B21__B11 = subtractMatrix(B21, B11, n2, p-p2);
    
    Matrix P4 = strassenRecursive(A22, B21__B11);
    
    auto A11_A12 = addMatrix(A11, A12, m-m2, n2);

    Matrix P5 = strassenRecursive(A11_A12, B22);
    
    auto A21__A11 = subtractMatrix(A21, A11, m-m2, n-n2);
    auto B11_B12 = addMatrix(B11, B12, n-n2, p-p2);
    
    Matrix P6 = strassenRecursive(A21__A11, B11_B12);

    auto A12__A22 = subtractMatrix(A12, A22, m-m2, n2);
    auto B21_B22 = addMatrix(B21, B22, n2, p-p2);

    Matrix P7 = strassenRecursive(A12__A22, B21_B22);

    // Składanie wynikowych bloków
    Matrix C11 = 
    addMatrix(
        subtractMatrix(
                addMatrix(P1, P4, m-m2, p-p2),
        P5, m-m2, p-p2),

    P7, m-m2, p-p2);



    Matrix C12 = addMatrix(P3, P5, m-m2, p2);
    Matrix C21 = addMatrix(P2, P4, m2, p-p2);

    Matrix C22 =
    addMatrix(
        subtractMatrix(
            addMatrix(P1, P3, m2, p2),
        P2, m2, p2),
    P6, m2, p2);



    // Wynikowa macierz
    Matrix C = createMatrix(m, p);
    paste(C, C11, 0, 0);
    paste(C, C12, 0, p-p2);
    paste(C, C21, m-m2, 0);
    paste(C, C22, m-m2, p-p2);

    return C;
}


// Funkcja pomocnicza — uruchamia mnożenie rekurencyjne
Matrix multiplyStrassen(const Matrix& A, const Matrix& B) {
    // Ma być tylko dla kwadratowych (choć w sumie chyba by działało bez tego)
    if (B.size() != A[0].size() || A[0].size() != A.size() || B.size() != B[0].size())
        throw std::invalid_argument("Niepoprawne rozmiary macierzy, mnożenie niemożliwe");

    
    return strassenRecursive(A, B);
}