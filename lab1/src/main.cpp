#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>

#include "multiply_binet.hpp"
#include "multiply_strassen.hpp"

#include <windows.h>
#include <psapi.h>

using namespace std;


// Tworzy losową macierz kwadratową
Matrix randomMatrix(int n, int seed = 0) {
    mt19937 gen(seed);
    uniform_real_distribution<double> dist(-5.0, 5.0);
    Matrix M(n, vector<double>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            M[i][j] = dist(gen);
    return M;
}

// Pomiar czasu w sekundach
template <typename Func>
double measureTime(Func f) {
    auto start = chrono::high_resolution_clock::now();
    f();
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end - start;
    return diff.count();
}

// Pomiar aktualnego zużycia pamięci (w MB)
double getMemoryUsageMB() {
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return static_cast<double>(info.WorkingSetSize) / (1024.0 * 1024.0);
}


int main() {
    // Zapis wyników do pliku
    ofstream out("results.csv");
    out << "n,time_seconds,mem_delta_MB\n";

    for (int n = 100; n <= 3000; n += 100) {
        cout << "Measuring n = " << n << "...\n";

        Matrix A = randomMatrix(n, n);
        Matrix B = randomMatrix(n, n);

        // Pomiar czasu i pamięci
        double mem_before = getMemoryUsageMB();

        double t1 = measureTime([&]() {
            auto C = multiplyBinet(A, B);
        });

        double mem_after = getMemoryUsageMB();
        double mem_delta = mem_after - mem_before;

        out << n << "," << fixed << setprecision(8) << t1
            << "," << setprecision(4) << mem_delta << "\n";
    }

    out.close();
    cout << "\nWyniki zapisano do pliku results.csv\n";
    return 0;
}
