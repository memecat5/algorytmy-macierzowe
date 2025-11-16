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
#endif
#include <vector>

using namespace std;

// Struktura do przechowywania macierzy i licznika operacji
struct Matrix {
    vector<vector<double>> data;
    int n;
    
    Matrix(int size) : n(size), data(size, vector<double>(size, 0.0)) {}
    
    void print(const string& name = "") const {
        if (!name.empty()) cout << name << ":\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
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

// 1. GENEROWANIE LOSOWYCH MACIERZY
Matrix generateRandomMatrix(int n) {
    Matrix m(n);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(1e-8, 1.0);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            m.data[i][j] = dis(gen);
        }
    }
    
    return m;
}

// Funkcja pomocnicza - kopiowanie fragmentu macierzy
void copySubmatrix(const Matrix& src, Matrix& dst, 
                   int srcRow, int srcCol, int dstRow, int dstCol, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            dst.data[dstRow + i][dstCol + j] = src.data[srcRow + i][srcCol + j];
        }
    }
}

// MNOŻENIE MACIERZY - KLASYCZNE O(n^3)
Matrix multiplyMatrices(const Matrix& A, const Matrix& B, OpCounter& counter) {
    int n = A.n;
    Matrix C(n);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C.data[i][j] = 0.0;
            for (int k = 0; k < n; k++) {
                C.data[i][j] += A.data[i][k] * B.data[k][j];
                counter.multiplications++;
                if (k > 0) counter.additions++;
            }
        }
    }
    
    return C;
}

// MNOŻENIE MACIERZY - REKURENCYJNE O(n^3)
Matrix multiplyRecursive(const Matrix& A, const Matrix& B, OpCounter& counter, int threshold = 32);

// Dodawanie macierzy
Matrix addMatrices(const Matrix& A, const Matrix& B, OpCounter& counter) {
    int n = A.n;
    Matrix C(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C.data[i][j] = A.data[i][j] + B.data[i][j];
            counter.additions++;
        }
    }
    return C;
}

// Odejmowanie macierzy
Matrix subtractMatrices(const Matrix& A, const Matrix& B, OpCounter& counter) {
    int n = A.n;
    Matrix C(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C.data[i][j] = A.data[i][j] - B.data[i][j];
            counter.additions++;
        }
    }
    return C;
}

Matrix multiplyRecursive(const Matrix& A, const Matrix& B, OpCounter& counter, int threshold) {
    int n = A.n;
    
    // Przypadek bazowy - użyj klasycznego algorytmu dla małych macierzy
    if (n <= threshold) {
        return multiplyMatrices(A, B, counter);
    }
    
    // Podział macierzy na ćwiartki
    int half = n / 2;
    
    Matrix A11(half), A12(half), A21(half), A22(half);
    Matrix B11(half), B12(half), B21(half), B22(half);
    
    // Kopiowanie bloków
    copySubmatrix(A, A11, 0, 0, 0, 0, half);
    copySubmatrix(A, A12, 0, half, 0, 0, half);
    copySubmatrix(A, A21, half, 0, 0, 0, half);
    copySubmatrix(A, A22, half, half, 0, 0, half);
    
    copySubmatrix(B, B11, 0, 0, 0, 0, half);
    copySubmatrix(B, B12, 0, half, 0, 0, half);
    copySubmatrix(B, B21, half, 0, 0, 0, half);
    copySubmatrix(B, B22, half, half, 0, 0, half);
    
    // Rekurencyjne mnożenie: 8 wywołań
    Matrix C11_1 = multiplyRecursive(A11, B11, counter, threshold);
    Matrix C11_2 = multiplyRecursive(A12, B21, counter, threshold);
    Matrix C11 = addMatrices(C11_1, C11_2, counter);
    
    Matrix C12_1 = multiplyRecursive(A11, B12, counter, threshold);
    Matrix C12_2 = multiplyRecursive(A12, B22, counter, threshold);
    Matrix C12 = addMatrices(C12_1, C12_2, counter);
    
    Matrix C21_1 = multiplyRecursive(A21, B11, counter, threshold);
    Matrix C21_2 = multiplyRecursive(A22, B21, counter, threshold);
    Matrix C21 = addMatrices(C21_1, C21_2, counter);
    
    Matrix C22_1 = multiplyRecursive(A21, B12, counter, threshold);
    Matrix C22_2 = multiplyRecursive(A22, B22, counter, threshold);
    Matrix C22 = addMatrices(C22_1, C22_2, counter);
    
    // Składanie wyniku
    Matrix C(n);
    copySubmatrix(C11, C, 0, 0, 0, 0, half);
    copySubmatrix(C12, C, 0, 0, 0, half, half);
    copySubmatrix(C21, C, 0, 0, half, 0, half);
    copySubmatrix(C22, C, 0, 0, half, half, half);
    
    return C;
}

// 2. REKURENCYJNE ODWRACANIE MACIERZY
Matrix invertMatrixRecursive(const Matrix& m, OpCounter& counter) {
    int n = m.n;
    
    if (n == 1) {
        Matrix inv(1);
        inv.data[0][0] = 1.0 / m.data[0][0];
        counter.divisions++;
        return inv;
    }
    
    int half = n / 2;
    int remainder = n - half;
    
    Matrix A(half), B(half), C(remainder), D(remainder);
    
    copySubmatrix(m, A, 0, 0, 0, 0, half);
    copySubmatrix(m, B, 0, half, 0, 0, half);
    copySubmatrix(m, C, half, 0, 0, 0, remainder);
    copySubmatrix(m, D, half, half, 0, 0, remainder);
    
    Matrix A_inv = invertMatrixRecursive(A, counter);
    Matrix CA_inv = multiplyMatrices(C, A_inv, counter);
    Matrix CA_invB = multiplyMatrices(CA_inv, B, counter);
    Matrix S = subtractMatrices(D, CA_invB, counter);
    Matrix S_inv = invertMatrixRecursive(S, counter);
    
    Matrix result(n);
    Matrix A_invB = multiplyMatrices(A_inv, B, counter);
    Matrix S_invC = multiplyMatrices(S_inv, C, counter);
    Matrix A_invB_S_inv = multiplyMatrices(A_invB, S_inv, counter);
    Matrix S_invC_A_inv = multiplyMatrices(S_invC, A_inv, counter);
    Matrix temp = multiplyMatrices(A_invB_S_inv, S_invC_A_inv, counter);
    
    for (int i = 0; i < half; i++) {
        for (int j = 0; j < half; j++) {
            result.data[i][j] = A_inv.data[i][j] + temp.data[i][j];
            counter.additions++;
        }
    }
    
    Matrix minusA_invB_S_inv(half);
    for (int i = 0; i < half; i++) {
        for (int j = 0; j < remainder; j++) {
            minusA_invB_S_inv.data[i][j] = -A_invB_S_inv.data[i][j];
        }
    }
    copySubmatrix(minusA_invB_S_inv, result, 0, 0, 0, half, half);
    
    Matrix minusS_invC_A_inv(remainder);
    for (int i = 0; i < remainder; i++) {
        for (int j = 0; j < half; j++) {
            minusS_invC_A_inv.data[i][j] = -S_invC_A_inv.data[i][j];
        }
    }
    copySubmatrix(minusS_invC_A_inv, result, 0, 0, half, 0, remainder);
    copySubmatrix(S_inv, result, 0, 0, half, half, remainder);
    
    return result;
}

// 3. REKURENCYJNA ELIMINACJA GAUSSA
void gaussEliminationRecursive(Matrix& m, Matrix& b, int start, int end, OpCounter& counter) {
    if (start >= end) return;
    
    int n = end - start;
    if (n == 1) return;
    
    for (int i = start + 1; i <= end; i++) {
        if (fabs(m.data[start][start]) < 1e-12) {
            throw runtime_error("Dzielenie przez zero w eliminacji Gaussa");
        }
        
        double factor = m.data[i][start] / m.data[start][start];
        counter.divisions++;
        
        for (int j = start; j <= end; j++) {
            m.data[i][j] -= factor * m.data[start][j];
            counter.multiplications++;
            counter.additions++;
        }
        
        for (int j = 0; j < b.n; j++) {
            b.data[i][j] -= factor * b.data[start][j];
            counter.multiplications++;
            counter.additions++;
        }
    }
    
    gaussEliminationRecursive(m, b, start + 1, end, counter);
}

Matrix solveGaussRecursive(Matrix m, Matrix b, OpCounter& counter) {
    int n = m.n;
    gaussEliminationRecursive(m, b, 0, n - 1, counter);
    
    Matrix x(n);
    for (int i = n - 1; i >= 0; i--) {
        for (int j = 0; j < b.n; j++) {
            x.data[i][j] = b.data[i][j];
            
            for (int k = i + 1; k < n; k++) {
                x.data[i][j] -= m.data[i][k] * x.data[k][j];
                counter.multiplications++;
                counter.additions++;
            }
            
            x.data[i][j] /= m.data[i][i];
            counter.divisions++;
        }
    }
    
    return x;
}

// 4. REKURENCYJNA FAKTORYZACJA LU
void luDecompositionRecursive(Matrix& m, int start, int end, OpCounter& counter) {
    if (start >= end) return;
    
    int n = end - start + 1;
    if (n == 1) return;
    
    for (int i = start + 1; i <= end; i++) {
        if (fabs(m.data[start][start]) < 1e-12) {
            throw runtime_error("Dzielenie przez zero w faktoryzacji LU");
        }
        m.data[i][start] /= m.data[start][start];
        counter.divisions++;
    }
    
    for (int i = start + 1; i <= end; i++) {
        for (int j = start + 1; j <= end; j++) {
            m.data[i][j] -= m.data[i][start] * m.data[start][j];
            counter.multiplications++;
            counter.additions++;
        }
    }
    
    luDecompositionRecursive(m, start + 1, end, counter);
}

double computeDeterminantFromLU(const Matrix& lu, OpCounter& counter) {
    double det = 1.0;
    for (int i = 0; i < lu.n; i++) {
        det *= lu.data[i][i];
        counter.multiplications++;
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
    Matrix identity = multiplyMatrices(m, inv, verifyCounter);
    
    double error = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            error += fabs(identity.data[i][j] - expected);
        }
    }
    cout << "Błąd weryfikacji: " << error << "\n";
}

void testGaussElimination(int n) {
    cout << "\n======================================\n";
    cout << "TEST ELIMINACJI GAUSSA (n=" << n << ")\n";
    cout << "======================================\n";
    
    Matrix A = generateRandomMatrix(n);
    Matrix b(n);
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
        luDecompositionRecursive(lu, 0, n - 1, counter);
        double det = computeDeterminantFromLU(lu, counter);
        
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
            luDecompositionRecursive(lu, 0, n - 1, counter4);
            computeDeterminantFromLU(lu, counter4);
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
    return 0;
#endif
}

// Benchmark: test multiplication for n = 1..maxN
void benchmarkAllSizes(int maxN = 1000) {
    cout << "\n============================================\n";
    cout << " BENCHMARK: MATRIX MULTIPLICATION (1.." << maxN << ")\n";
    cout << "============================================\n";

    ofstream file("matrix_benchmark.csv");
    file << "n,time_ms,FLOPs,memory_MB,peak_memory_MB,FLOPs_per_sec\n";

    for (int n = 1; n <= maxN; n++) {
        Matrix A = generateRandomMatrix(n);
        Matrix B = generateRandomMatrix(n);

        // Theoretical FLOPs for classic O(n³) multiplication
        double flops = 2.0 * pow(n, 3); // n³ multiplications + n³ additions

        // Estimate memory used by A, B, and result C (3 * n² doubles)
        double memoryUsageMB = (3.0 * n * n * sizeof(double)) / (1024.0 * 1024.0);

        OpCounter counter;
        auto mem_before = getMemoryUsageMB();
        auto start = chrono::high_resolution_clock::now();

        Matrix C = multiplyMatrices(A, B, counter);

        auto end = chrono::high_resolution_clock::now();
        auto mem_after = getMemoryUsageMB();
        double time_ms = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;

        double peakMem = max(mem_before, mem_after);
        double flopsPerSec = (flops / (time_ms / 1000.0));

        file << n << ","
             << time_ms << ","
             << flops << ","
             << memoryUsageMB << ","
             << peakMem << ","
             << flopsPerSec << "\n";

        if (n % 100 == 0 || n == 1 || n == maxN) {
            cout << setw(5) << n
                 << " | time=" << setw(8) << fixed << setprecision(4) << time_ms << " ms"
                 << " | mem=" << setw(7) << setprecision(3) << memoryUsageMB << " MB"
                 << " | FLOPs/s=" << scientific << flopsPerSec << fixed << "\n";
        }
    }

    file.close();
    cout << "\nBenchmark completed. Results saved to 'matrix_benchmark.csv'\n";
}

void showMenu() {
    cout << "\n╔════════════════════════════════════════════╗\n";
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