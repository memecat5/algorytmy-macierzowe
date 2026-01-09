#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>

// Dołączenie biblioteki Eigen
#include <Eigen/Dense>
#include <Eigen/SVD>

// Implementacja bibliotek stb
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace Eigen;
using namespace std;

// --- STRUKTURY DANYCH ---

struct HNode {
    bool is_leaf;
    int t_min, t_max, s_min, s_max;
    int rank;
    
    MatrixXd U;
    VectorXd SingularValues;
    MatrixXd V; 

    std::vector<HNode*> children;

    HNode(int tmin, int tmax, int smin, int smax) 
        : t_min(tmin), t_max(tmax), s_min(smin), s_max(smax), is_leaf(false), rank(0) {}

    ~HNode() {
        for (auto child : children) delete child;
    }
};

// --- FUNKCJE POMOCNICZE (GetBlock, SVD, Draw) ---

MatrixXd GetBlock(const MatrixXd& full_matrix, int t_min, int t_max, int s_min, int s_max) {
    return full_matrix.block(t_min, s_min, t_max - t_min + 1, s_max - s_min + 1);
}

// Rekurencyjna kompresja (zgodnie ze slajdem 32 i 25)
HNode* CompressMatrix(const MatrixXd& A, int t_min, int t_max, int s_min, int s_max, int max_rank, double epsilon) {
    HNode* node = new HNode(t_min, t_max, s_min, s_max);
    int rows = t_max - t_min + 1;
    int cols = s_max - s_min + 1;

    // Warunek stopu dla pojedynczych pikseli
    if (rows <= 1 || cols <= 1) {
        node->is_leaf = true;
        JacobiSVD<MatrixXd> svd(GetBlock(A, t_min, t_max, s_min, s_max), ComputeThinU | ComputeThinV);
        node->rank = svd.rank();
        node->U = svd.matrixU();
        node->SingularValues = svd.singularValues();
        node->V = svd.matrixV();
        return node;
    }

    MatrixXd block = GetBlock(A, t_min, t_max, s_min, s_max);
    
    // Używamy BDCSVD dla większych bloków (szybsze), Jacobi dla małych
    BDCSVD<MatrixXd> svd(block, ComputeThinU | ComputeThinV);
    const VectorXd& sing_vals = svd.singularValues();
    
    // Warunek dopuszczalności (Admissibility Condition)
    // Sprawdzamy ile wartości osobliwych jest znaczących (> epsilon)
    int k = 0;
    for (int i = 0; i < sing_vals.size(); ++i) {
        if (sing_vals(i) > epsilon) k++;
    }

    // Heurystyka: jeśli ranga jest mała (<= max_rank) ORAZ mniejsza niż połowa boku (slajd 26)
    bool admissible = (k <= max_rank) && (k <= std::min(rows, cols) / 2);

    if (admissible) {
        node->is_leaf = true;
        node->rank = k;
        if (k > 0) {
            node->U = svd.matrixU().leftCols(k);
            node->SingularValues = sing_vals.head(k);
            node->V = svd.matrixV().leftCols(k);
        }
    } else {
        // Podział na 4 synów
        int t_mid = t_min + (rows / 2) - 1;
        int s_mid = s_min + (cols / 2) - 1;
        if (t_mid < t_min) t_mid = t_min;
        if (s_mid < s_min) s_mid = s_min;

        node->children.push_back(CompressMatrix(A, t_min, t_mid, s_min, s_mid, max_rank, epsilon));
        if (s_mid + 1 <= s_max)
            node->children.push_back(CompressMatrix(A, t_min, t_mid, s_mid + 1, s_max, max_rank, epsilon));
        if (t_mid + 1 <= t_max)
            node->children.push_back(CompressMatrix(A, t_mid + 1, t_max, s_min, s_mid, max_rank, epsilon));
        if (t_mid + 1 <= t_max && s_mid + 1 <= s_max)
            node->children.push_back(CompressMatrix(A, t_mid + 1, t_max, s_mid + 1, s_max, max_rank, epsilon));
    }
    return node;
}

void DecompressMatrix(HNode* node, MatrixXd& Output) {
    if (node->is_leaf) {
        if (node->rank > 0) {
            MatrixXd Block = node->U * node->SingularValues.asDiagonal() * node->V.transpose();
            Output.block(node->t_min, node->s_min, 
                         node->t_max - node->t_min + 1, 
                         node->s_max - node->s_min + 1) = Block;
        } else {
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

void DrawStructure(HNode* node, MatrixXd& GridMap, double val = 0.0) {
    if (node->is_leaf) {
        int r1 = node->t_min; int r2 = node->t_max;
        int c1 = node->s_min; int c2 = node->s_max;
        // Rysuj ramkę
        for (int c = c1; c <= c2; ++c) { GridMap(r1, c) = val; GridMap(r2, c) = val; }
        for (int r = r1; r <= r2; ++r) { GridMap(r, c1) = val; GridMap(r, c2) = val; }
    } else {
        for (auto child : node->children) {
            DrawStructure(child, GridMap, val);
        }
    }
}

// --- FUNKCJA RAPORTUJĄCA ---
// Uruchamia jeden przypadek testowy, zapisuje obraz i strukturę
void RunTestCase(const std::string& caseName, 
                 const MatrixXd& R, const MatrixXd& G, const MatrixXd& B, 
                 int width, int height, 
                 int r_param, double epsilon) {
    
    std::cout << ">>> Uruchamianie: " << caseName 
              << " (r=" << r_param << ", epsilon=" << epsilon << ")" << std::endl;

    // 1. Kompresja
    HNode* rootR = CompressMatrix(R, 0, height - 1, 0, width - 1, r_param, epsilon);
    HNode* rootG = CompressMatrix(G, 0, height - 1, 0, width - 1, r_param, epsilon);
    HNode* rootB = CompressMatrix(B, 0, height - 1, 0, width - 1, r_param, epsilon);

    // 2. Dekompresja
    MatrixXd R_out(height, width), G_out(height, width), B_out(height, width);
    DecompressMatrix(rootR, R_out);
    DecompressMatrix(rootG, G_out);
    DecompressMatrix(rootB, B_out);

    // 3. Zapis wyniku (Obraz)
    std::vector<unsigned char> output_data(width * height * 3);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int idx = (i * width + j) * 3;
            output_data[idx]     = (unsigned char)std::min(255.0, std::max(0.0, R_out(i, j)));
            output_data[idx + 1] = (unsigned char)std::min(255.0, std::max(0.0, G_out(i, j)));
            output_data[idx + 2] = (unsigned char)std::min(255.0, std::max(0.0, B_out(i, j)));
        }
    }
    std::string filenameImg = "out_" + caseName + ".jpg";
    stbi_write_jpg(filenameImg.c_str(), width, height, 3, output_data.data(), 90);

    // 4. Zapis wyniku (Struktura - na bazie kanału R)
    MatrixXd Grid(height, width);
    Grid.fill(255.0); // Białe tło
    DrawStructure(rootR, Grid, 0.0); // Czarne linie
    
    std::vector<unsigned char> grid_data(width * height);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
             grid_data[i * width + j] = (unsigned char)Grid(i, j);
        }
    }
    std::string filenameGrid = "grid_" + caseName + ".jpg";
    stbi_write_jpg(filenameGrid.c_str(), width, height, 1, grid_data.data(), 90);

    // Czyszczenie
    delete rootR; delete rootG; delete rootB;
    std::cout << "    Zapisano: " << filenameImg << " oraz " << filenameGrid << std::endl;
}

int main() {
    const char* filename = "input.jpg"; 
    int width, height, channels;
    unsigned char* img_data = stbi_load(filename, &width, &height, &channels, 3);

    if (!img_data) {
        std::cerr << "Blad: Nie mozna otworzyc input.jpg" << std::endl;
        return -1;
    }
    std::cout << "Wczytano obraz: " << width << "x" << height << std::endl;

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

    // --- ZADANIE 3: Wykres wartości osobliwych (Singular Values) ---
    // Slajd 33: "Proszę uruchomić SVD dla całej bitmapy R... znaleźć wartości osobliwe"
    std::cout << "Obliczanie pelnego SVD dla kanalu R (moze potrwac)..." << std::endl;
    BDCSVD<MatrixXd> full_svd(R, ComputeThinU | ComputeThinV);
    VectorXd sigmas = full_svd.singularValues();

    // Zapis wartości do CSV (do stworzenia wykresu w Excelu/Pythonie)
    std::ofstream csvFile("singular_values.csv");
    for (int i = 0; i < sigmas.size(); ++i) {
        csvFile << i << "," << sigmas(i) << "\n";
    }
    csvFile.close();
    std::cout << "Zapisano wartosci osobliwe do singular_values.csv" << std::endl;

    // Pobranie kluczowych wartości osobliwych do metod kompresji
    double sigma_1 = sigmas(0);                   // Największa
    double sigma_N = sigmas(sigmas.size() - 1);   // Najmniejsza (sigma_2^k)
    double sigma_mid = sigmas(sigmas.size() / 2); // Środkowa (sigma_2^k / 2)

    std::cout << "Sigma_1 (max): " << sigma_1 << std::endl;
    std::cout << "Sigma_mid: " << sigma_mid << std::endl;
    std::cout << "Sigma_N (min): " << sigma_N << std::endl;

    // --- IMPLEMENTACJA METOD ZE SLAJDU 33 ---
    
    // Metoda 1: r=1, delta = sqrt(sigma_1) 
    RunTestCase("Metoda1_r1_sqrtSigma1", R, G, B, width, height, 1, std::sqrt(sigma_1));

    // Metoda 2: r=1, delta = sigma_N 
    // (Bardzo dokładna, mały próg błędu)
    RunTestCase("Metoda2_r1_SigmaN", R, G, B, width, height, 1, sigma_N);

    // Metoda 3: r=1, delta = sigma_mid 
    RunTestCase("Metoda3_r1_SigmaMid", R, G, B, width, height, 1, sigma_mid);

    // Metoda 4: r=4, delta = sigma_1 
    // (Bardzo duży próg błędu -> duża kompresja, niska jakość)
    RunTestCase("Metoda4_r4_Sigma1", R, G, B, width, height, 4, sigma_1);

    // Metoda 5: r=4, delta = sigma_N 
    RunTestCase("Metoda5_r4_SigmaN", R, G, B, width, height, 4, sigma_N);

    // Metoda 6: r=4, delta = sigma_mid 
    RunTestCase("Metoda6_r4_SigmaMid", R, G, B, width, height, 4, sigma_mid);

    // Metoda Dodatkowa (Twoja własna, optymalna) [cite: 538]
    // Spróbujmy dobrać parametry "rozsądne" wizualnie
    double optimal_eps = sigma_1 * 0.02; // np. 2% największej wartości
    RunTestCase("MetodaOptymalna", R, G, B, width, height, 8, optimal_eps);

    std::cout << "Wszystkie testy zakonczone." << std::endl;
    return 0;
}