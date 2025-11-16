#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>

using namespace std;

// Struktura do przechowywania macierzy i licznika operacji
struct Matrix {
    vector<vector<double>> data;
    int rows, cols;
    
    Matrix(int n) : rows(n), cols(n), data(n, vector<double>(n, 0.0)) {}
    Matrix(int rows, int cols) : rows(rows), cols(cols), data(rows, vector<double>(cols, 0.0)) {}
    
    // Operator dodawania macierzy
    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::runtime_error("Niezgodne wymiary macierzy (dodawanie)");
        }
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return result;
    }

    // Operator odejmowania macierzy
    Matrix operator-(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::runtime_error("Niezgodne wymiary macierzy (odejmowanie)");
        }
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] - other.data[i][j];
            }
        }
        return result;
    }

    // Operator mnożenia przez skalar
    Matrix operator*(double scalar) const {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] * scalar;
            }
        }
        return result;
    }


    void print(const string& name = "") const {
        if (!name.empty()) cout << name << ":\n";
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                cout << setw(12) << fixed << setprecision(6) << data[i][j] << " ";
            }
            cout << "\n";
        }
        cout << "\n";
    }
};

// Licznik operacji zmiennoprzecinkowych
struct OpCounter {
    long long additions = 0;
    long long multiplications = 0;
    long long divisions = 0;
    
    void reset() {
        additions = multiplications = divisions = 0;
    }
    
    long long total() const {
        return additions + multiplications + divisions;
    }
    
    void print(const string& operation) const {
        cout << "=== Statystyki operacji dla: " << operation << " ===\n";
        cout << "Dodawania/odejmowania: " << additions << "\n";
        cout << "Mnożenia: " << multiplications << "\n";
        cout << "Dzielenia: " << divisions << "\n";
        cout << "Razem: " << total() << "\n\n";
    }
};

Matrix createMatrix(int rows, int cols, double value = 0.);
Matrix multiplyBinet(const Matrix& A, const Matrix& B, OpCounter& counter);

// Tworzy pustą macierz o wymiarach rows x cols
Matrix createMatrix(int rows, int cols, double value) {
    return Matrix(rows, cols);
}

// Rekurencyjne mnożenie fragmentów macierzy
void multiplyRecursive(
    const Matrix& A, int a_row, int a_col,
    const Matrix& B, int b_row, int b_col,
    Matrix& C, int c_row, int c_col,
    int m, int n, int p,
    OpCounter& counter)
{   
    // x_col, x_row - lewy górny róg rozważanej podmacierzy

    // m - liczba wierszy A, n - liczba kolumn A (i wierszy B), p - liczba kolumn B

    // przypadek brzegowy - mnożenie samych liczb
    if (m == 1 && n == 1 && p == 1) {
        C.data[c_row][c_col] += A.data[a_row][a_col] * B.data[b_row][b_col];
        ++counter.additions;
        ++counter.multiplications;
        return;
    }

    // Przy dzieleniu może wyjść wymiar równy zero - to znaczy że
    // dzieliliśmy macierz, której któryś z wymiarów był równy jeden.
    // To jest odpowiednik mnożenia przez 0 przy paddingu.
    if (m == 0 || n == 0 || p == 0)
        return;

    // jeśli któryś wymiar = 1 klasyczne mnożenie
    if (m <= 1 || n <= 1 || p <= 1) {
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < p; ++j)
                for (int k = 0; k < n; ++k){
                    C.data[c_row + i][c_col + j] += A.data[a_row + i][a_col + k] * B.data[b_row + k][b_col + j],
                    ++counter.additions;
                    ++counter.multiplications;
                }
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

    std::pair<long long, long long> recursive_count;
    // A11 * B11
    multiplyRecursive(A, a11_row, a11_col, B, b11_row, b11_col, C, c11_row, c11_col, m2, n2, p2, counter);
    // A12 * B21
    multiplyRecursive(A, a12_row, a12_col, B, b21_row, b21_col, C, c11_row, c11_col, m2, n-n2, p2, counter);


    // --------------------------------C12--------------------------------

    // A11 * B12
    multiplyRecursive(A, a11_row, a11_col, B, b12_row, b12_col, C, c12_row, c12_col, m2, n2, p-p2, counter);

    // A12 * B22
    multiplyRecursive(A, a12_row, a12_col, B, b22_row, b22_col, C, c12_row, c12_col, m2, n-n2, p-p2, counter);
    
    
    // --------------------------------C21--------------------------------


    // A21 * B11
    multiplyRecursive(A, a21_row, a21_col, B, b11_row, b11_col, C, c21_row, c21_col, m-m2, n2, p2, counter);
    
    // A22 * B21
    multiplyRecursive(A, a22_row, a22_col, B, b21_row, b21_col, C, c21_row, c21_col, m-m2, n-n2, p2, counter);

    // --------------------------------C22--------------------------------
    
    // A21 * B12
    multiplyRecursive(A, a21_row, a21_col, B, b12_row, b12_col, C, c22_row, c22_col, m-m2, n2, p-p2, counter); 

    // A22 * B22
    multiplyRecursive(A, a22_row, a22_col, B, b22_row, b22_col, C, c22_row, c22_col, m-m2, n-n2, p-p2, counter);

    return;
}

// Funkcja pomocnicza — uruchamia mnożenie rekurencyjne
Matrix multiplyBinet(const Matrix& A, const Matrix& B, OpCounter& counter) {
    int m = A.data.size();          // wiersze A
    int n = A.data[0].size();       // kolumny A = wiersze B
    int p = B.data[0].size();       // kolumny B

    if (B.data.size() != A.data[0].size())
        throw std::invalid_argument("Niepoprawne rozmiary macierzy, mnozenie niemozliwe");

    Matrix C = createMatrix(m, p);
    multiplyRecursive(A, 0, 0, B, 0, 0, C, 0, 0, m, n, p, counter);

    return C;
}