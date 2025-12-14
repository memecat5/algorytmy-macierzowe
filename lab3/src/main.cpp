#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <fstream>
#include <random>
#include <memory>

// Zakładamy, że pliki nagłówkowe są w tym samym katalogu
#include <Eigen/Core>
#include <Eigen/Dense>
#include "compress.hpp"   // Definicja SVDNode i compressMatrixRecursive
#include "partial_svd.hpp"
#include "bitmap_utils.hpp"

using namespace Eigen;
using namespace std;

void drawCompressedMatrix(const SVDNode* node, Ref<MatrixXd> target) {
    if (node->isLeaf) {
        target = node->U * node->S.asDiagonal() * node->V.transpose();
    } else {
        int r = target.rows(), c = target.cols();
        drawCompressedMatrix(node->children[0].get(), target.block(0, 0, r/2, c/2));
        drawCompressedMatrix(node->children[1].get(), target.block(0, c/2, r/2, c - c/2));
        drawCompressedMatrix(node->children[2].get(), target.block(r/2, 0, r - r/2, c/2));
        drawCompressedMatrix(node->children[3].get(), target.block(r/2, c/2, r - r/2, c - c/2));
    }
}

int interactiveDrawCompressedMatrix() {
    int W = 512, H = 512;
    MatrixXd R, G, B;

    // wygenerowanie bitmapy (ale takiej słabej, nie polecam)
    /*
    generateBitmap(W, H, R, G, B);
    savePPM("original.ppm", R, G, B);
    */

    // Wczytanie bitmapy z pliku
    loadPPM("original.ppm", R, G, B);

    int b ;//= 20;        
    double delta;// = 50000.0; 

    cout << "Podaj parametry dokladnosci - delta: ";
    cin >> delta;
    cout << "b: ";
    cin >> b;

    cout << "Start (Spectra SVD 1.0+, b=" << b << ")..." << endl;

    vector<std::unique_ptr<SVDNode>> channels(3);
    channels[0] = compressMatrixRecursive(R, delta, b);
    channels[1] = compressMatrixRecursive(G, delta, b);
    channels[2] = compressMatrixRecursive(B, delta, b);

    MatrixXd oR(H,W), oG(H,W), oB(H,W);
    drawCompressedMatrix(channels[0].get(), oR);
    drawCompressedMatrix(channels[1].get(), oG);
    drawCompressedMatrix(channels[2].get(), oB);

    savePPM("final.ppm", oR, oG, oB);
    cout << "Zapisano final.ppm" << endl;
    return 0;
}

// Funkcja generująca macierz o topologii siatki 3D zgodnie z poleceniem
// "wiersz = wierzchołek, niezerowe losowe wartości w kolumnach sąsiadujące wierzchołki siatki"
MatrixXd generate3DGridMatrix(int k) {
    int n = std::pow(2, k); // bok sześcianu
    int N = n * n * n;      // całkowity rozmiar macierzy (N x N)
    
    // Używamy macierzy gęstej, bo algorytm kompresji w compress.cpp przyjmuje MatrixXd
    MatrixXd A = MatrixXd::Zero(N, N);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.1, 1.0); // Losowe wartości

    // Iteracja po wszystkich punktach siatki (z, y, x)
    for (int z = 0; z < n; ++z) {
        for (int y = 0; y < n; ++y) {
            for (int x = 0; x < n; ++x) {
                // Indeks wiersza odpowiadający punktowi (x,y,z)
                // Mapowanie [i,j,k] -> row (zgodnie ze slajdem 8/27, wzór podobny)
                int row = z * n * n + y * n + x;
                
                // Sąsiedzi w 3D: 6 kierunków
                int dx[] = {1, -1, 0, 0, 0, 0};
                int dy[] = {0, 0, 1, -1, 0, 0};
                int dz[] = {0, 0, 0, 0, 1, -1};

                for (int i = 0; i < 6; ++i) {
                    int nx = x + dx[i];
                    int ny = y + dy[i];
                    int nz = z + dz[i];

                    // Sprawdzenie czy sąsiad jest wewnątrz siatki
                    if (nx >= 0 && nx < n && ny >= 0 && ny < n && nz >= 0 && nz < n) {
                        int col = nz * n * n + ny * n + nx;
                        A(row, col) = dis(gen); 
                    }
                }
                // Wartość na przekątnej (opcjonalnie, dla stabilności numerycznej)
                A(row, row) = 1.0; 
            }
        }
    }
    return A;
}

// Rekurencyjne mnożenie macierzy skompresowanej przez wektor
// Zgodnie z algorytmem na slajdzie 20/27
VectorXd h_matrix_vector_mult(const SVDNode* node, const VectorXd& x) {
    if (node->isLeaf) {
        if (node->S.size() == 0) return VectorXd::Zero(x.rows());
        return node->U * (node->S.asDiagonal() * (node->V.transpose() * x));
    } else {
        int half = x.rows() / 2;
        VectorXd x1 = x.head(half);
        VectorXd x2 = x.tail(x.rows() - half);
        VectorXd y1 = h_matrix_vector_mult(node->children[0].get(), x1) + 
                      h_matrix_vector_mult(node->children[1].get(), x2);
        VectorXd y2 = h_matrix_vector_mult(node->children[2].get(), x1) + 
                      h_matrix_vector_mult(node->children[3].get(), x2);
        VectorXd y(x.rows());
        y << y1, y2;
        return y;
    }
}

int main() {
    ofstream results("results.csv");
    // Dodano kolumnę SSE_Error (błąd mnożenia) oraz Matrix_Diff_Norm (błąd rekonstrukcji macierzy)
    results << "k,N,Time_Compression,Time_MVM,MVM_FLOPs,MMM_FLOPs,SSE_Error,Matrix_Diff_Norm\n";
    
    double delta = 1e-6;
    int b = 10;

    for (int k : {2, 3, 4}) {
        int n_side = std::pow(2, k);
        int N = n_side * n_side * n_side;
        
        cout << "--- Przetwarzanie k=" << k << " (N=" << N << "x" << N << ") ---" << endl;

        cout << "Generowanie macierzy..." << endl;
        MatrixXd A = generate3DGridMatrix(k);
        
        cout << "Kompresja..." << endl;
        auto start = chrono::high_resolution_clock::now();
        auto root = compressMatrixRecursive(A, delta, b);
        auto end = chrono::high_resolution_clock::now();
        double time_comp = chrono::duration<double>(end - start).count();
        cout << "Czas kompresji: " << time_comp << " s" << endl;

        // --- Analiza złożoności ---
        long long mvm_flops = getMVMFlops(root.get());
        long long mmm_flops = getMMMFlops(root.get(), root.get());

        // --- Mnożenie Macierz-Wektor ---
        VectorXd x = VectorXd::Random(N); 
        
        // 1. Dokładne mnożenie (Dense)
        VectorXd y_dense = A * x;

        // 2. Mnożenie H-Macierzy (Approx)
        cout << "Mnożenie macierz-wektor..." << endl;
        start = chrono::high_resolution_clock::now();
        VectorXd y_approx = h_matrix_vector_mult(root.get(), x);
        end = chrono::high_resolution_clock::now();
        double time_mvm = chrono::duration<double>(end - start).count();

        // 3. Obliczenie błędu SSE (Sum of Squared Errors)
        // SSE = ||y_dense - y_approx||^2
        double sse = (y_dense - y_approx).squaredNorm();
        cout << "Błąd MVM (SSE): " << sse << endl;

        // --- Rekonstrukcja macierzy i sprawdzenie normy błędu ---
        // (Opcjonalne, ale potwierdza, że dekonstrukcja działa)
        MatrixXd A_rec = decompressMatrix(root.get());
        double matrix_diff_norm = (A - A_rec).norm(); // Norma Frobeniusa różnicy
        cout << "Błąd rekonstrukcji macierzy (Norma Frobeniusa): " << matrix_diff_norm << endl;

        // Zapis do CSV
        results << k << "," << N << "," << time_comp << "," << time_mvm << "," 
                << mvm_flops << "," << mmm_flops << "," << sse << "," << matrix_diff_norm << "\n";
        results.flush();
    }

    results.close();
    cout << "Zakończono." << endl;
    return 0;
}