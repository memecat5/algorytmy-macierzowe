#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#ifndef _WIN32
    #include <sys/resource.h>  // for memory usage (Linux/macOS)
    #include <unistd.h>
#else
    #include <windows.h>
    #include <psapi.h>
#endif
#include <vector>
#include "binet.hpp"

using namespace std;

// Helper: build random matrix with fixed seed for reproducibility
Matrix generateRandomMatrix(int n) {
    std::mt19937 rng(314159);
    std::uniform_real_distribution<double> dist(-2.0, 2.0);

    Matrix M(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            M.data[i][j] = dist(rng);
    return M;
}

void splitMatrix(const Matrix& A, Matrix& A11, Matrix& A12, Matrix& A21, Matrix& A22) {
    int k = A11.rows; // lub A11.cols
    int m = A22.rows; // lub A22.cols

    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            A11.data[i][j] = A.data[i][j];
        }
        for (int j = 0; j < m; ++j) {
            A12.data[i][j] = A.data[i][j + k];
        }
    }
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < k; ++j) {
            A21.data[i][j] = A.data[i + k][j];
        }
        for (int j = 0; j < m; ++j) {
            A22.data[i][j] = A.data[i + k][j + k];
        }
    }
}

/**
 * Funkcja pomocnicza: Łączy 4 bloki w B
 * B11 (k x k), B12 (k x m), B21 (m x k), B22 (m x m)
 */
void joinMatrix(Matrix& B, const Matrix& B11, const Matrix& B12, const Matrix& B21, const Matrix& B22) {
    int k = B11.rows; // lub B11.cols
    int m = B22.rows; // lub B22.cols

    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            B.data[i][j] = B11.data[i][j];
        }
        for (int j = 0; j < m; ++j) {
            B.data[i][j + k] = B12.data[i][j];
        }
    }
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < k; ++j) {
            B.data[i + k][j] = B21.data[i][j];
        }
        for (int j = 0; j < m; ++j) {
            B.data[i + k][j + k] = B22.data[i][j];
        }
    }
}

Matrix invertMatrixRecursive(const Matrix& A, OpCounter& counter) {
    if (A.rows != A.cols) {
        throw std::runtime_error("Macierz nie jest kwadratowa.");
    }
    if (A.rows == 0) return Matrix(0, 0);

    // BAZA REKURENCJI
    if (A.rows == 1) {
        Matrix B(1, 1);
        if (std::fabs(A.data[0][0]) < 1e-12) {
            throw std::runtime_error("Dzielenie przez zero: Macierz osobliwa.");
        }
        B.data[0][0] = 1.0 / A.data[0][0];
        ++counter.divisions;
        return B;
    }
    
    // 1. PODZIAŁ
    int n = A.rows;
    int k = n / 2;    // Rozmiar A11 (podłoga)
    int m = n - k;    // Rozmiar A22 (sufit)

    Matrix A11(k, k), A12(k, m), A21(m, k), A22(m, m);
    splitMatrix(A, A11, A12, A21, A22);

    // 2. REKURENCJA
    
    // Oblicz A11_inv = A11^{-1}
    Matrix A11_inv = invertMatrixRecursive(A11, counter);

    // Oblicz S = A22 - A21 * A11_inv * A12
    Matrix S_tmp = multiplyBinet(A21, A11_inv, counter);
    Matrix S = A22 - multiplyBinet(S_tmp, A12, counter);
    counter.additions += A22.rows * A22.cols;

    // Oblicz S_inv = S^{-1}
    Matrix S_inv = invertMatrixRecursive(S, counter);

    // Teraz mamy A11_inv i S_inv. Obliczamy bloki B.
    // B22 = S_inv
    Matrix B22 = S_inv;

    // B12 = -A11_inv * A12 * B22
    Matrix B12_tmp = multiplyBinet(A11_inv, A12, counter);
    Matrix B12 = multiplyBinet(B12_tmp, B22, counter) * -1.0;
    counter.multiplications += B12.rows * B12.cols;

    // B21 = -B22 * A21 * A11_inv
    Matrix B21_tmp = multiplyBinet(B22, A21, counter);
    Matrix B21 = multiplyBinet(B21_tmp, A11_inv, counter) * -1.0;
    counter.multiplications += B21.rows * B12.cols;

    // B11 = A11_inv - (A11_inv * A12 * B21)
    Matrix B11_tmp = multiplyBinet(A11_inv, A12, counter);
    Matrix B11 = A11_inv - multiplyBinet(B11_tmp, B21, counter);
    counter.additions -= A11_inv.rows * A11_inv.cols;

    // 3. POŁĄCZENIE
    Matrix B(n, n);
    joinMatrix(B, B11, B12, B21, B22);

    return B;
}

/**
 * 4. NOWE: Funkcje pomocnicze do dzielenia/łączenia WEKTORÓW
 * (Wektory traktujemy jako macierze n x 1)
 */
void splitVector(const Matrix& v, Matrix& v1, Matrix& v2) {
    int k = v1.rows; // Rozmiar v1
    int m = v2.rows; // Rozmiar v2
    for (int i = 0; i < k; ++i) {
        v1.data[i][0] = v.data[i][0];
    }
    for (int i = 0; i < m; ++i) {
        v2.data[i][0] = v.data[i + k][0];
    }
}

void joinVector(Matrix& v, const Matrix& v1, const Matrix& v2) {
    int k = v1.rows;
    int m = v2.rows;
    for (int i = 0; i < k; ++i) {
        v.data[i][0] = v1.data[i][0];
    }
    for (int i = 0; i < m; ++i) {
        v.data[i + k][0] = v2.data[i][0];
    }
}


/**
 * 5. NOWA GŁÓWNA FUNKCJA: Rekurencyjny blokowy solver Gaussa
 * Rozwiązuje A*x = b
 */
Matrix solveGaussRecursive(const Matrix& A, const Matrix& b, OpCounter& counter) {
    if (A.rows != A.cols) throw std::runtime_error("Macierz A nie jest kwadratowa.");
    if (A.rows != b.rows || b.cols != 1) throw std::runtime_error("Niezgodne wymiary A i b.");

    int n = A.rows;
    if (n == 0) return Matrix(0, 1);

    // BAZA REKURENCJI (n = 1)
    if (n == 1) {
        if (std::fabs(A.data[0][0]) < 1e-12) throw std::runtime_error("Dzielenie przez zero.");
        Matrix x(1, 1);
        x.data[0][0] = b.data[0][0] / A.data[0][0];
        ++counter.divisions;
        return x;
    }

    // 1. PODZIAŁ
    int k = n / 2;
    int m = n - k;

    // Podział A
    Matrix A11(k, k), A12(k, m), A21(m, k), A22(m, m);
    splitMatrix(A, A11, A12, A21, A22);

    // Podział b
    Matrix b1(k, 1), b2(m, 1);
    splitVector(b, b1, b2);

    // 2. OBLICZENIA (Kroki eliminacji)
    
    // Potrzebujemy A11_inv do obliczenia S i b2'
    Matrix A11_inv = invertMatrixRecursive(A11, counter);

    // Oblicz S = A22 - A21 * A11_inv * A12
    Matrix S = A22 - multiplyBinet(multiplyBinet(A21, A11_inv, counter), A12, counter);
    counter.additions += A22.cols * A22.rows;

    // Oblicz b2' = b2 - A21 * A11_inv * b1
    Matrix b2_prime = b2 - multiplyBinet(multiplyBinet(A21, A11_inv, counter), b1, counter);
    counter.additions += b2.rows * b2.cols;

    // 3. REKURENCJA

    // Krok 1: Rozwiąż S * x2 = b2'
    Matrix x2 = solveGaussRecursive(S, b2_prime, counter);

    // Krok 2: Rozwiąż A11 * x1 = b1' (gdzie b1' = b1 - A12 * x2)
    Matrix b1_prime = b1 - multiplyBinet(A12, x2, counter);
    counter.additions += b1.rows * b1.cols;

    Matrix x1 = solveGaussRecursive(A11, b1_prime, counter);

    // 4. POŁĄCZENIE
    Matrix x(n, 1);
    joinVector(x, x1, x2);
    
    return x;
}

/**
 * Funkcja pomocnicza: Kopiuje blok 'source' do 'target' z offsetem
 */
void copyBlock(Matrix& target, const Matrix& source, int row_offset, int col_offset) {
    for (int i = 0; i < source.rows; ++i) {
        for (int j = 0; j < source.cols; ++j) {
            target.data[i + row_offset][j + col_offset] = source.data[i][j];
        }
    }
}

/**
 * 1. NOWY HELPER: Rozwiązuje LX = B (dla kroku U12)
 * L jest dolnotrójkątna (z 1 na diagonali), B ma k kolumn
 */
Matrix solveL(const Matrix& L, const Matrix& B, OpCounter& counter) {
    int n = L.rows;
    int k = B.cols;
    Matrix X(n, k);

    for (int j = 0; j < k; ++j) { // Dla każdej kolumny B
        for (int i = 0; i < n; ++i) { // Podstawienie w przód
            double sum = 0.0;
            for (int p = 0; p < i; ++p) {
                sum += L.data[i][p] * X.data[p][j];
                ++counter.additions;
            }
            // L.data[i][i] to 1.0 (zgodnie z Doolittle)
            X.data[i][j] = (B.data[i][j] - sum);
            ++counter.additions;
        }
    }
    return X;
}

/**
 * 2. NOWY HELPER: Rozwiązuje XU = B (dla kroku L21)
 * U jest górnotrójkątna, B ma k wierszy
 */
Matrix solveU_T(const Matrix& U, const Matrix& B, OpCounter& counter) {
    int n = U.cols;
    int k = B.rows;
    Matrix X(k, n);

    for (int i = 0; i < k; ++i) { // Dla każdego wiersza B
        for (int j = n - 1; j >= 0; --j) { // Podstawienie wsteczne (wierszowe)
            double sum = 0.0;
            for (int p = j + 1; p < n; ++p) {
                sum += X.data[i][p] * U.data[p][j];
                ++counter.additions;
            }
            if (std::fabs(U.data[j][j]) < 1e-12) {
                throw std::runtime_error("Dzielenie przez zero w solveU_T");
            }
            X.data[i][j] = (B.data[i][j] - sum) / U.data[j][j];
            ++counter.additions;
            ++counter.divisions;
        }
    }
    return X;
}

/**
 * GŁÓWNA FUNKCJA: Rekurencyjna blokowa faktoryzacja LU
 * Zwraca parę {L, U}
 */
std::pair<Matrix, Matrix> blockLU(const Matrix& A, OpCounter& counter) {
    if (A.rows != A.cols) {
        throw std::runtime_error("Macierz nie jest kwadratowa.");
    }

    int n = A.rows;
    Matrix L(n, n);
    Matrix U(n, n);

    // BAZA REKURENCJI (n = 1) - Rozkład Doolittle'a
    if (n == 1) {
        if (std::fabs(A.data[0][0]) < 1e-12) {
             // Zwykle tutaj robimy pivoting, ale w tym algorytmie...
            throw std::runtime_error("Element [0][0] jest zerowy. Wymagany pivoting.");
        }
        L.data[0][0] = 1.0;
        U.data[0][0] = A.data[0][0];
        return {L, U};
    }

    // 1. PODZIAŁ
    int k = n / 2;
    int m = n - k;
    Matrix A11(k, k), A12(k, m), A21(m, k), A22(m, m);
    splitMatrix(A, A11, A12, A21, A22);

    // 2. REKURENCJA
    
    // Krok 1: A11 = L11 * U11
    auto [L11, U11] = blockLU(A11, counter);

    // Krok 2: U12 = L11^{-1} * A12
    Matrix U12 = solveL(L11, A12, counter);

    // Krok 3: L21 = A21 * U11^{-1}
    Matrix L21 = solveU_T(U11, A21, counter);

    // Krok 4: S = A22 - L21 * U12
    Matrix S = A22 - multiplyBinet(L21, U12, counter);

    // Krok 5: S = L22 * U22
    auto [L22, U22] = blockLU(S, counter);

    // 3. POŁĄCZENIE
    // Złożenie macierzy L
    copyBlock(L, L11, 0, 0);       // L[0:k, 0:k] = L11
    copyBlock(L, L21, k, 0);       // L[k:n, 0:k] = L21
    copyBlock(L, L22, k, k);       // L[k:n, k:n] = L22
    // Blok L12 (górny prawy) pozostaje zerowy

    // Złożenie macierzy U
    copyBlock(U, U11, 0, 0);       // U[0:k, 0:k] = U11
    copyBlock(U, U12, 0, k);       // U[0:k, k:n] = U12
    copyBlock(U, U22, k, k);       // U[k:n, k:n] = U22
    // Blok U21 (dolny lewy) pozostaje zerowy

    return {L, U};
}

double detLU(const Matrix& A, OpCounter& counter){
    auto [L, U] = blockLU(A, counter);
    double det = 1.;
    for(int i = 0; i < U.rows; ++i){
        det *= U.data[i][i];
    }
    return det;
}

void testInversion(int n) {
    cout << "\n======================================\n";
    cout << "TEST ODWRACANIA MACIERZY (n=" << n << ")\n";
    cout << "======================================\n";
    
    Matrix m = generateRandomMatrix(n);
    if (n <= 5) m.print("Macierz wejściowa");
    
    OpCounter counter;
    Matrix inv = invertMatrixRecursive(m, counter);
    
    if (n <= 5) inv.print("Macierz odwrotna");
    counter.print("Odwracanie macierzy");
    
    OpCounter verifyCounter;
    Matrix identity = multiplyBinet(m, inv, verifyCounter);
    
    double error = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            error += fabs(identity.data[i][j] - expected);
        }
    }
    cout << "Blad weryfikacji: " << error << "\n";
}

void testGaussElimination(int n) {
    cout << "\n======================================\n";
    cout << "TEST ELIMINACJI GAUSSA (n=" << n << ")\n";
    cout << "======================================\n";
    
    Matrix A = generateRandomMatrix(n);
    Matrix b(n, 1);
    for (int i = 0; i < n; i++) {
        b.data[i][0] = 1.0;
    }
    
    if (n <= 5) {
        A.print("Macierz A");
        b.print("Wektor b");
    }
    
    OpCounter counter;
    Matrix x = solveGaussRecursive(A, b, counter);
    
    if (n <= 5) x.print("Rozwiązanie x");
    counter.print("Eliminacja Gaussa");
}

void testLUFactorization(int n) {
    cout << "\n======================================\n";
    cout << "TEST FAKTORYZACJI LU (n=" << n << ")\n";
    cout << "======================================\n";
    
    Matrix m = generateRandomMatrix(n);
    if (n <= 5) m.print("Macierz wejściowa");
    
    Matrix lu = m;
    OpCounter counter;
    
    try {
        blockLU(lu, counter);
        double det = detLU(lu, counter);
        
        if (n <= 5) lu.print("Macierz LU");
        cout << "Wyznacznik: " << scientific << det << fixed << "\n";
        counter.print("Faktoryzacja LU + wyznacznik");
    } catch (const exception& e) {
        cout << "Błąd: " << e.what() << "\n";
    }
}

// Wykres dla wszystkich operacji
void plotAllOperations(const vector<int>& sizes) {
    vector<double> x_axis;
    vector<double> invert_ops, gauss_ops, lu_ops;

    cout << "\nGenerowanie danych dla wszystkich operacji...\n";

    for (int n : sizes) {
        cout << "Testowanie dla n = " << n << "..." << endl;
        x_axis.push_back(n);

        Matrix A = generateRandomMatrix(n);
        Matrix B = generateRandomMatrix(n);

        // Odwracanie
        OpCounter counter2;
        try {
            invertMatrixRecursive(A, counter2);
            invert_ops.push_back(counter2.total());
        } catch (...) {
            invert_ops.push_back(0);
        }

        // Gauss
        Matrix b(n);
        for (int i = 0; i < n; i++) b.data[i][0] = 1.0;
        OpCounter counter3;
        try {
            solveGaussRecursive(A, b, counter3);
            gauss_ops.push_back(counter3.total());
        } catch (...) {
            gauss_ops.push_back(0);
        }

        // LU
        Matrix lu = A;
        OpCounter counter4;
        try {
            blockLU(lu, counter4);
            detLU(lu, counter4);
            lu_ops.push_back(counter4.total());
        } catch (...) {
            lu_ops.push_back(0);
        }
    }

    // Zapis do CSV
    ofstream file("all_operations_comparison.csv");
    file << "n,multiply_ops,invert_ops,gauss_ops,lu_ops\n";
    for (size_t i = 0; i < x_axis.size(); ++i) {
        file << x_axis[i] << ","
             << invert_ops[i] << ","
             << gauss_ops[i] << ","
             << lu_ops[i] << "\n";
    }
    file.close();

    cout << "\nDane zapisane do pliku 'all_operations_comparison.csv'\n";
}

// Helper: get current memory usage in MB
double getMemoryUsageMB() {
#ifndef _WIN32
        struct rusage usage{};
        getrusage(RUSAGE_SELF, &usage);
    #ifdef __APPLE__
        return usage.ru_maxrss / (1024.0 * 1024.0);
    #else
        return usage.ru_maxrss / 1024.0;  // Linux returns kB
    #endif
#else
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return static_cast<double>(info.WorkingSetSize) / (1024.0 * 1024.0);
#endif
}

void benchmarkInversion(int maxN, int step){
    ofstream file("inversion_benchmark.csv");
    file << "n,time_ms,FLOPs,peak_memory_MB,FLOPs_per_sec\n";

    for (int n = 1; n <= maxN; n += step) {
        Matrix A = generateRandomMatrix(n);

        OpCounter counter;
        auto mem_before = getMemoryUsageMB();
        auto start = chrono::high_resolution_clock::now();

        Matrix C = invertMatrixRecursive(A, counter);

        auto end = chrono::high_resolution_clock::now();
        auto mem_after = getMemoryUsageMB();
        double time_ms = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;

        double peakMem = max(mem_before, mem_after);    // chyba tylko mem after?
        double flopsPerSec = (counter.total() / (time_ms / 1000.0));

        file << n << ","
             << time_ms << ","
             << counter.total() << ","
             << peakMem << ","
             << flopsPerSec << "\n";
    }

    file.close();
}

void benchmarkGauss(int maxN, int step){
    ofstream file("gauss_benchmark.csv");
    file << "n,time_ms,FLOPs,peak_memory_MB,FLOPs_per_sec\n";

    for (int n = 1; n <= maxN; n += step) {
        Matrix A = generateRandomMatrix(n);
        Matrix B(n, 1);

        for(int i=0; i<n; ++i){
            B.data[i][0] = 1.;
        }

        OpCounter counter;
        auto mem_before = getMemoryUsageMB();
        auto start = chrono::high_resolution_clock::now();

        Matrix C = solveGaussRecursive(A, B, counter);

        auto end = chrono::high_resolution_clock::now();
        auto mem_after = getMemoryUsageMB();
        double time_ms = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;

        double peakMem = max(mem_before, mem_after);    // chyba tylko mem after?
        double flopsPerSec = (counter.total() / (time_ms / 1000.0));

        file << n << ","
             << time_ms << ","
             << counter.total() << ","
             << peakMem << ","
             << flopsPerSec << "\n";
    }

    file.close();
}

void benchmarkLUdecomposition(int maxN, int step){
    ofstream file("lu_benchmark.csv");
    file << "n,time_ms,FLOPs,peak_memory_MB,FLOPs_per_sec\n";

    for (int n = 1; n <= maxN; n += step) {
        Matrix A = generateRandomMatrix(n);

        OpCounter counter;
        auto mem_before = getMemoryUsageMB();
        auto start = chrono::high_resolution_clock::now();

        double det = detLU(A, counter);

        auto end = chrono::high_resolution_clock::now();
        auto mem_after = getMemoryUsageMB();
        double time_ms = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;

        double peakMem = max(mem_before, mem_after);    // chyba tylko mem after?
        double flopsPerSec = (counter.total() / (time_ms / 1000.0));

        file << n << ","
             << time_ms << ","
             << counter.total() << ","
             << peakMem << ","
             << flopsPerSec << "\n";
    }

    file.close();
}

// Benchmark: test all operations for n = 1..maxN
void benchmarkAllSizes(int maxN = 1000, int step = 10) {
    cout << "\n"
         << "============================================\n"
         << " BENCHMARK: (1.." << maxN << ")\n"
         << "============================================\n";

    benchmarkInversion(maxN, step);
    benchmarkGauss(maxN, step);
    benchmarkLUdecomposition(maxN, step);

    cout << "Benchmark completed, results saved to corresponding csv files\n";

}


void showMenu() {
    cout << "\n"
         << "╔════════════════════════════════════════════╗\n";
    cout << "║   OPERACJE NA MACIERZACH - MENU GŁÓWNE     ║\n";
    cout << "╚════════════════════════════════════════════╝\n";
    cout << "1. Porównanie algorytmów mnożenia macierzy\n";
    cout << "2. Odwracanie macierzy\n";
    cout << "3. Eliminacja Gaussa\n";
    cout << "4. Faktoryzacja LU i wyznacznik\n";
    cout << "5. Testy dla wszystkich operacji\n";
    cout << "6. Wykresy zależności operacji od rozmiaru macierzy\n";
    cout << "7. Benchmark (1..1000)\n";
    cout << "0. Wyjście\n";
    cout << "\nWybierz opcję: ";
}

// Funkcja do generowania wykresów
void generatePlots() {
    cout << "\n╔════════════════════════════════════════════╗\n";
    cout << "║        GENEROWANIE WYKRESÓW                ║\n";
    cout << "╚════════════════════════════════════════════╝\n";
    cout << "Wybierz typ wykresu:\n";
    cout << "1. Porównanie algorytmów mnożenia\n";
    cout << "2. Wszystkie operacje macierzowe\n";
    cout << "3. Szczegółowa analiza mnożenia (czas vs operacje)\n";
    cout << "\nWybór: ";
    
    int plotChoice;
    cin >> plotChoice;
    
    vector<int> sizes;
    cout << "\nPodaj rozmiary macierzy do testowania (potęgi 2, np: 4 8 16 32 64)\n";
    cout << "Zakończ wpisując 0: ";
    
    int size;
    while (cin >> size && size > 0) {
        sizes.push_back(size);
    }
    
    if (sizes.empty()) {
        cout << "Brak podanych rozmiarów!\n";
        return;
    }
    
    sort(sizes.begin(), sizes.end());
    
    switch (plotChoice) {
        case 1:
            //plotMultiplicationComparison(sizes);
            cout << "Option removed\n";
            break;
        case 2:
            plotAllOperations(sizes);
            break;
        case 3:
            //plotDetailedMultiplication(sizes);
            cout << "Option removed\n";
            break;
        default:
            cout << "Nieprawidłowy wybór!\n";
    }
}

int main() {
    cout << fixed << setprecision(6);
    
    int choice;
    int n;
    
    while (true) {
        showMenu();
        cin >> choice;
        
        if (choice == 0) {
            cout << "Dziękuję za użycie programu!\n";
            break;
        }
        
        if (choice >= 1 && choice <= 4) {
            cout << "Podaj rozmiar macierzy (potęga 2 zalecana, np. 4, 8, 16, 32): ";
            cin >> n;
            
            if (n <= 0 || n > 1024) {
                cout << "Nieprawidłowy rozmiar macierzy!\n";
                continue;
            }
        }
        
        switch (choice) {
            case 1:
                //testMultiplication(n);
                cout << "Mutiplication removed\n";
                break;
            case 2:
                testInversion(n);
                break;
            case 3:
                testGaussElimination(n);
                break;
            case 4:
                testLUFactorization(n);
                break;
            case 5:
                cout << "Podaj rozmiar macierzy: ";
                cin >> n;
                if (n <= 0 || n > 1024) {
                    cout << "Nieprawidłowy rozmiar macierzy!\n";
                    continue;
                }
                //testMultiplication(n);
                cout << "Mutiplication removed\n";
                testInversion(n);
                testGaussElimination(n);
                testLUFactorization(n);
                break;
            case 6:
                generatePlots();
                break;
            case 7:
                benchmarkAllSizes();
                break;
            default:
                cout << "Nieprawidłowa opcja!\n";
        }
    }
    
    return 0;
}