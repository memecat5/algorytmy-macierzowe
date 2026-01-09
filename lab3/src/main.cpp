#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

// Dołączenie biblioteki Eigen
#include <Eigen/Dense>
#include <Eigen/SVD>

// Implementacja bibliotek stb (odkomentuj jeśli masz pliki nagłówkowe)
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace Eigen;
using namespace std;

// Struktura węzła drzewa (H-Matrix Node)
// Bazuje na strukturze ze slajdu 25 [cite: 417]
struct HNode {
    bool is_leaf;
    int t_min, t_max, s_min, s_max; // Zakresy indeksów
    int rank;
    
    // Dane dla liścia (SVD)
    MatrixXd U;
    VectorXd SingularValues;
    MatrixXd V; // Przechowujemy V, rekonstrukcja to U * S * V.transpose()

    // Dzieci (4 ćwiartki)
    std::vector<HNode*> children;

    HNode(int tmin, int tmax, int smin, int smax) 
        : t_min(tmin), t_max(tmax), s_min(smin), s_max(smax), is_leaf(false), rank(0) {}

    ~HNode() {
        for (auto child : children) delete child;
    }
};

// Funkcja pomocnicza do pobrania bloku macierzy
MatrixXd GetBlock(const MatrixXd& full_matrix, int t_min, int t_max, int s_min, int s_max) {
    return full_matrix.block(t_min, s_min, t_max - t_min + 1, s_max - s_min + 1);
}

// Główna funkcja rekurencyjna: CreateTree / CompressMatrix
// Implementacja algorytmu ze slajdu 32 [cite: 494, 495]
HNode* CompressMatrix(const MatrixXd& A, int t_min, int t_max, int s_min, int s_max, int max_rank, double epsilon) {
    HNode* node = new HNode(t_min, t_max, s_min, s_max);
    int rows = t_max - t_min + 1;
    int cols = s_max - s_min + 1;

    // Warunek stopu dla bardzo małych bloków (np. 1 piksel)
    if (rows <= 1 || cols <= 1) {
        node->is_leaf = true;
        // Pełne SVD dla małego bloku
        JacobiSVD<MatrixXd> svd(GetBlock(A, t_min, t_max, s_min, s_max), ComputeThinU | ComputeThinV);
        node->rank = svd.rank();
        node->U = svd.matrixU();
        node->SingularValues = svd.singularValues();
        node->V = svd.matrixV();
        return node;
    }

    // 1. Wykonaj SVD dla bloku (używamy BDCSVD dla większych macierzy, Jacobi dla mniejszych)
    // Slajd 32: [U, D, V] = truncatedSVD(block, r+1) 
    MatrixXd block = GetBlock(A, t_min, t_max, s_min, s_max);
    BDCSVD<MatrixXd> svd(block, ComputeThinU | ComputeThinV);
    
    const VectorXd& sing_vals = svd.singularValues();
    
    // Sprawdzenie warunku dopuszczalności (Admissibility Condition)
    // Slajd 32 pkt 2: if D(r+1, r+1) < epsilon 
    // Musimy sprawdzić, czy odrzucając wartości osobliwe powyżej max_rank, błąd jest akceptowalny (< epsilon)
    // ORAZ czy w ogóle opłaca się kompresować (heurystyka ze slajdu 26: k <= size/2) [cite: 427]
    
    int k = 0;
    for (int i = 0; i < sing_vals.size(); ++i) {
        if (sing_vals(i) > epsilon) k++;
    }

    bool admissible = (k <= max_rank) && (k <= std::min(rows, cols) / 2);

    if (admissible) {
        // KOMPRESJA (LEAF)
        node->is_leaf = true;
        node->rank = k;
        
        // Zapisujemy tylko k najważniejszych wartości (Truncated SVD)
        // Slajd 25 pkt 8-10 
        if (k > 0) {
            node->U = svd.matrixU().leftCols(k);
            node->SingularValues = sing_vals.head(k);
            node->V = svd.matrixV().leftCols(k);
        } else {
            // Blok zerowy lub pomijalny
            node->rank = 0;
        }
    } else {
        // PODZIAŁ (REKURENCJA)
        // Slajd 32 pkt 9 [cite: 504, 506]
        int t_mid = t_min + (rows / 2) - 1;
        int s_mid = s_min + (cols / 2) - 1;

        // Upewnij się, że podział jest poprawny
        if (t_mid < t_min) t_mid = t_min;
        if (s_mid < s_min) s_mid = s_min;

        // Tworzenie 4 synów (ćwiartki)
        // Lewy-Górny
        node->children.push_back(CompressMatrix(A, t_min, t_mid, s_min, s_mid, max_rank, epsilon));
        // Prawy-Górny
        if (s_mid + 1 <= s_max)
            node->children.push_back(CompressMatrix(A, t_min, t_mid, s_mid + 1, s_max, max_rank, epsilon));
        // Lewy-Dolny
        if (t_mid + 1 <= t_max)
            node->children.push_back(CompressMatrix(A, t_mid + 1, t_max, s_min, s_mid, max_rank, epsilon));
        // Prawy-Dolny
        if (t_mid + 1 <= t_max && s_mid + 1 <= s_max)
            node->children.push_back(CompressMatrix(A, t_mid + 1, t_max, s_mid + 1, s_max, max_rank, epsilon));
    }

    return node;
}

// Funkcja dekompresji - odtwarza macierz z drzewa
void DecompressMatrix(HNode* node, MatrixXd& Output) {
    if (node->is_leaf) {
        if (node->rank > 0) {
            // M = U * S * V^T
            MatrixXd Block = node->U * node->SingularValues.asDiagonal() * node->V.transpose();
            
            // Wstawienie bloku w odpowiednie miejsce
            Output.block(node->t_min, node->s_min, 
                         node->t_max - node->t_min + 1, 
                         node->s_max - node->s_min + 1) = Block;
        } else {
            // Blok zerowy
            Output.block(node->t_min, node->s_min, 
                         node->t_max - node->t_min + 1, 
                         node->s_max - node->s_min + 1).setZero();
        }
    } else {
        for (auto child : node->children) {
            DecompressMatrix(child, Output);
        }
    }
}

// "Rysowacz" struktury kompresji (siatka) - Zadanie 3 
// Rysuje ramki liści na białym tle
void DrawStructure(HNode* node, MatrixXd& GridMap, double val = 0.0) {
    if (node->is_leaf) {
        // Rysuj ramkę wokół bloku
        int r1 = node->t_min;
        int r2 = node->t_max;
        int c1 = node->s_min;
        int c2 = node->s_max;

        // Poziome linie
        for (int c = c1; c <= c2; ++c) {
            GridMap(r1, c) = val;
            GridMap(r2, c) = val;
        }
        // Pionowe linie
        for (int r = r1; r <= r2; ++r) {
            GridMap(r, c1) = val;
            GridMap(r, c2) = val;
        }
    } else {
        for (auto child : node->children) {
            DrawStructure(child, GridMap, val);
        }
    }
}

int main() {
    // 1. Wczytanie obrazu
    const char* filename = "../lena.jpg"; // Podmień na swój plik
    int width, height, channels;
    unsigned char* img_data = stbi_load(filename, &width, &height, &channels, 3); // Wymuszamy 3 kanały (RGB)

    if (!img_data) {
        std::cerr << "Nie udalo sie wczytac obrazu!" << std::endl;
        return -1;
    }

    std::cout << "Wczytano obraz: " << width << "x" << height << std::endl;

    // 2. Konwersja Bitmapy na 3 Macierze (R, G, B) [cite: 487]
    MatrixXd R(height, width), G(height, width), B(height, width);

    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int idx = (i * width + j) * 3;
            R(i, j) = (double)img_data[idx];
            G(i, j) = (double)img_data[idx + 1];
            B(i, j) = (double)img_data[idx + 2];
        }
    }
    stbi_image_free(img_data);

    // Obliczanie globalnych wartości osobliwych dla celów raportu [cite: 519]
    // Uwaga: Dla dużych obrazów pełne SVD jest kosztowne. Tutaj robimy to dla testu.
    BDCSVD<MatrixXd> full_svd(R, ComputeThinU | ComputeThinV); // SVD tylko dla kanału R jako przykład
    VectorXd singular_vals = full_svd.singularValues();
    double max_sv = singular_vals(0);
    std::cout << "Max Singular Value (Channel R): " << max_sv << std::endl;


    // 3. Ustawienie parametrów kompresji [cite: 490]
    // Przykładowe parametry (eksperyment 4 ze slajdu 33 [cite: 532])
    int r_param = 4;        // Max rank (b)
    double epsilon = max_sv * 0.05; // Próg odcięcia (np. 5% max wartości)

    std::cout << "Kompresja z parametrami: r=" << r_param << ", eps=" << epsilon << "..." << std::endl;

    // 4. Rekurencyjna kompresja (dla każdego kanału) [cite: 489]
    HNode* rootR = CompressMatrix(R, 0, height - 1, 0, width - 1, r_param, epsilon);
    HNode* rootG = CompressMatrix(G, 0, height - 1, 0, width - 1, r_param, epsilon);
    HNode* rootB = CompressMatrix(B, 0, height - 1, 0, width - 1, r_param, epsilon);

    std::cout << "Kompresja zakonczona." << std::endl;

    // 5. Rekonstrukcja (dekompresja) macierzy
    MatrixXd R_out(height, width), G_out(height, width), B_out(height, width);
    DecompressMatrix(rootR, R_out);
    DecompressMatrix(rootG, G_out);
    DecompressMatrix(rootB, B_out);

    // 6. Zapis skompresowanej bitmapy ("Rysowacz bitmapy") 
    std::vector<unsigned char> output_data(width * height * 3);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int idx = (i * width + j) * 3;
            // Clamp values to [0, 255]
            output_data[idx]     = (unsigned char)std::min(255.0, std::max(0.0, R_out(i, j)));
            output_data[idx + 1] = (unsigned char)std::min(255.0, std::max(0.0, G_out(i, j)));
            output_data[idx + 2] = (unsigned char)std::min(255.0, std::max(0.0, B_out(i, j)));
        }
    }

    stbi_write_jpg("output_compressed.jpg", width, height, 3, output_data.data(), 90);
    std::cout << "Zapisano output_compressed.jpg" << std::endl;

    // 7. Rysowanie struktury kompresji ("Rysowacz macierzy") 
    // Tworzymy białe tło
    MatrixXd Grid(height, width);
    Grid.fill(255.0);
    
    // Rysujemy czarne ramki (0.0) na podstawie drzewa kanału R
    DrawStructure(rootR, Grid, 0.0);

    std::vector<unsigned char> grid_data(width * height);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
             grid_data[i * width + j] = (unsigned char)Grid(i, j);
        }
    }
    stbi_write_jpg("../output_structure.jpg", width, height, 1, grid_data.data(), 90);
    std::cout << "Zapisano output_structure.jpg" << std::endl;

    // Czyszczenie pamięci
    delete rootR;
    delete rootG;
    delete rootB;

    return 0;
}