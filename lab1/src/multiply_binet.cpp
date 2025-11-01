#include <vector>

using Matrix = std::vector<std::vector<double>>;

bool equalMatrix(const Matrix& A, const Matrix& B) {
    if (A.size() != B.size() || A[0].size() != B[0].size()) return false;
    for (size_t i = 0; i < A.size(); ++i)
        for (size_t j = 0; j < A[0].size(); ++j)
            if (A[i][j] != B[i][j])
                return false;
    return true;
}

// Tworzy pustą macierz o wymiarach rows x cols
Matrix createMatrix(int rows, int cols, double value = 0) {
    return Matrix(rows, std::vector<double>(cols, value));
}

// Rekurencyjne mnożenie fragmentów macierzy
void multiplyRecursive(
    const Matrix& A, int a_row, int a_col,
    const Matrix& B, int b_row, int b_col,
    Matrix& C, int c_row, int c_col,
    int m, int n, int p)
{

    // x_col, x_row - lewy górny róg rozważanej podmacierzy

    // m - liczba wierszy A, n - liczba kolumn A (i wierszy B), p - liczba kolumn B

    // przypadek brzegowy - mnożenie samych liczb
    if (m == 1 && n == 1 && p == 1) {
        C[c_row][c_col] += A[a_row][a_col] * B[b_row][b_col];
        return;
    }

    // Przy dzieleniu może wyjść wymiar równy zero - to znaczy że
    // dzieliliśmy macierz, której któryś z wymiarów był równy jeden.
    // To jest odpowiednik mnożenia przez 0 przy paddingu.
    if (m == 0 || n == 0 || p == 0)
        return;

    // jeśli któryś wymiar = 1, możesz zamiast schodzić dalej zrobić klasyczne mnożenie
    // (ja w sumie nie rozumiem czemu to jest potrzebne, ale bez tego testy się wywalają, a czat mówi że ma być XD)
    if (m <= 1 || n <= 1 || p <= 1) {
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < p; ++j)
                for (int k = 0; k < n; ++k)
                    C[c_row + i][c_col + j] += A[a_row + i][a_col + k] * B[b_row + k][b_col + j];
        return;
    }
    
    
    // Dzielimy każdy wymiar na pół (dzielenie intów, więc tu może wyjść 0)
    int m2 = m / 2;
    int n2 = n / 2;
    int p2 = p / 2;

    /*         
               n2   n-n2
               ^     ^
        m2{   A11 | A12
    A =       ---------
        m-m2{ A21 | A22
    
               p2   p-p2
               ^     ^
        n2{   B11 | B12
    B =       ---------
        n-n2{ B21 | B22
    */

    // Dla czytelności, kompilator i tak to wywali
    
    // Podział A
    int a11_row = a_row;
    int a11_col = a_col;

    int a12_row = a_row;
    int a12_col = a_col + n2;

    int a21_row = a_row + m2;
    int a21_col = a_col;

    int a22_row = a_row + m2;
    int a22_col = a_col + n2;

    // Podział B
    int b11_row = b_row;
    int b11_col = b_col;

    int b12_row = b_row;
    int b12_col = b_col + p2;

    int b21_row = b_row + n2;
    int b21_col = b_col;

    int b22_row = b_row + n2;
    int b22_col = b_col + p2;

    // Podział C
    int c11_row = c_row;
    int c11_col = c_col;

    int c12_row = c_row;
    int c12_col = c_col + p2;

    int c21_row = c_row + m2;
    int c21_col = c_col;

    int c22_row = c_row + m2;
    int c22_col = c_col + p2;
    
    // --------------------------------C11--------------------------------

    // A11 * B11
    multiplyRecursive(A, a11_row, a11_col, B, b11_row, b11_col, C, c11_row, c11_col, m2, n2, p2);
    // A12 * B21
    multiplyRecursive(A, a12_row, a12_col, B, b21_row, b21_col, C, c11_row, c11_col, m2, n-n2, p2);


    // --------------------------------C12--------------------------------

    // A11 * B12
    multiplyRecursive(A, a11_row, a11_col, B, b12_row, b12_col, C, c12_row, c12_col, m2, n2, p-p2);
    // A12 * B22
    multiplyRecursive(A, a12_row, a12_col, B, b22_row, b22_col, C, c12_row, c12_col, m2, n-n2, p-p2);
    
    
    // --------------------------------C21--------------------------------


    // A21 * B11
    multiplyRecursive(A, a21_row, a21_col, B, b11_row, b11_col, C, c21_row, c21_col, m-m2, n2, p2);
    
    // A22 * B21
    multiplyRecursive(A, a22_row, a22_col, B, b21_row, b21_col, C, c21_row, c21_col, m-m2, n-n2, p2);

    // --------------------------------C22--------------------------------
    
    // A21 * B12
    multiplyRecursive(A, a21_row, a21_col, B, b12_row, b12_col, C, c22_row, c22_col, m-m2, n2, p-p2); 
    
    // A22 * B22
    multiplyRecursive(A, a22_row, a22_col, B, b22_row, b22_col, C, c22_row, c22_col, m-m2, n-n2, p-p2);

}

// Funkcja pomocnicza — uruchamia mnożenie rekurencyjne
Matrix multiplyMatrix(const Matrix& A, const Matrix& B) {
    int m = A.size();          // wiersze A
    int n = A[0].size();       // kolumny A = wiersze B
    int p = B[0].size();       // kolumny B

    Matrix C = createMatrix(m, p);
    multiplyRecursive(A, 0, 0, B, 0, 0, C, 0, 0, m, n, p);
    return C;
}