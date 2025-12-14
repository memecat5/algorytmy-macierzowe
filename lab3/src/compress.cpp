#include "compress.hpp"
#include "partial_svd.hpp"
#include <iostream>

std::unique_ptr<SVDNode> compressMatrixRecursive(const MatrixXd& block, double delta, int b) {
    auto node = std::make_unique<SVDNode>();
    int rows = block.rows();
    int cols = block.cols();

    // Próba obliczenia częściowego SVD
    MatrixXd U_part, V_part;
    VectorXd S_part;
    computePartialSVD(block, b, U_part, S_part, V_part);

    double smallestCalculatedSigma = (S_part.size() > 0) ? S_part(S_part.size() - 1) : 0.0;
    
    // Sprawdzenie warunku kompresji
    bool goodApproximation = (S_part.size() < b) || (smallestCalculatedSigma < delta);

    if (goodApproximation) {
        node->isLeaf = true;
        int keep = 0;
        for(int i = 0; i < S_part.size(); ++i) {
            if (S_part(i) >= delta) keep++;
        }
        if (keep == 0) keep = 1;

        node->U = U_part.leftCols(keep);
        node->S = S_part.head(keep);
        node->V = V_part.leftCols(keep);
    } else {
        node->isLeaf = false;
        int hRows = rows / 2;
        int hCols = cols / 2;
        // Kolejność: 0:TL, 1:TR, 2:BL, 3:BR
        node->children.push_back(compressMatrixRecursive(block.block(0, 0, hRows, hCols), delta, b));
        node->children.push_back(compressMatrixRecursive(block.block(0, hCols, hRows, cols - hCols), delta, b));
        node->children.push_back(compressMatrixRecursive(block.block(hRows, 0, rows - hRows, hCols), delta, b));
        node->children.push_back(compressMatrixRecursive(block.block(hRows, hCols, rows - hRows, cols - hCols), delta, b));
    }

    return node;
}

// --- Implementacja funkcji analitycznych ---

int getRows(const SVDNode* node) {
    if (node->isLeaf) return node->U.rows();
    // Suma wierszy górnych bloków (TL + BL to błąd, TL i TR mają tyle samo wierszy)
    // W układzie blokowym wiersze to suma wierszy (TL) + wierszy (BL)
    // children[0] to TL, children[2] to BL
    return getRows(node->children[0].get()) + getRows(node->children[2].get());
}

int getCols(const SVDNode* node) {
    if (node->isLeaf) return node->V.rows();
    // Kolumny to suma kolumn (TL) + kolumn (TR)
    // children[0] to TL, children[1] to TR
    return getCols(node->children[0].get()) + getCols(node->children[1].get());
}

long long getMVMFlops(const SVDNode* node) {
    if (!node) return 0;

    if (node->isLeaf) {
        long long k = node->S.size();
        long long m = node->U.rows(); 
        long long n = node->V.rows(); 

        // Koszt y = U * S * V^T * x
        // 1. V^T * x (n*1) -> k iloczynów skalarnych długości n -> 2*k*n
        // 2. Skalowanie S -> k mnożeń -> k
        // 3. U * wektor (k*1) -> m iloczynów skalarnych długości k -> 2*m*k
        return 2 * k * n + k + 2 * m * k;
    } else {
        long long flops = 0;
        for (const auto& child : node->children) {
            flops += getMVMFlops(child.get());
        }
        // Dodatkowo koszt scalenia wyników (dodawanie wektorów y1+y2)
        // Pomijamy dla uproszczenia, bo dominują mnożenia
        return flops;
    }
}

long long getMMMFlops(const SVDNode* A, const SVDNode* B) {
    if (!A || !B) return 0;

    // Przypadek 1: A jest liściem (Niska Ranga), B dowolne
    if (A->isLeaf) {
        // A = U_A * S_A * V_A^T. Mnożenie A * B = U_A * S_A * (V_A^T * B)
        // Kluczowy koszt to pomnożenie wektorów z V_A^T przez macierz B.
        // V_A^T ma wymiar k_A x colsA. 
        // Traktujemy to jak k_A operacji MVM na macierzy B.
        long long kA = A->S.size();
        return kA * getMVMFlops(B);
    }

    // Przypadek 2: A węzeł, B jest liściem (Niska Ranga)
    if (B->isLeaf) {
        // A * B = (A * U_B) * S_B * V_B^T
        // Kluczowy koszt to A * U_B. 
        // U_B ma wymiar rowsB x k_B.
        // Traktujemy to jak k_B operacji MVM na macierzy A.
        long long kB = B->S.size();
        return kB * getMVMFlops(A);
    }

    // Przypadek 3: A i B to węzły wewnętrzne -> Rekurencja
    // C = A * B
    // C00 = A00*B00 + A01*B10
    // C01 = A00*B01 + A01*B11
    // C10 = A10*B00 + A11*B10
    // C11 = A10*B01 + A11*B11
    
    long long totalFlops = 0;
    
    // 8 mnożeń rekurencyjnych
    totalFlops += getMMMFlops(A->children[0].get(), B->children[0].get()); // A00 * B00
    totalFlops += getMMMFlops(A->children[1].get(), B->children[2].get()); // A01 * B10
    
    totalFlops += getMMMFlops(A->children[0].get(), B->children[1].get()); // A00 * B01
    totalFlops += getMMMFlops(A->children[1].get(), B->children[3].get()); // A01 * B11
    
    totalFlops += getMMMFlops(A->children[2].get(), B->children[0].get()); // A10 * B00
    totalFlops += getMMMFlops(A->children[3].get(), B->children[2].get()); // A11 * B10
    
    totalFlops += getMMMFlops(A->children[2].get(), B->children[1].get()); // A10 * B01
    totalFlops += getMMMFlops(A->children[3].get(), B->children[3].get()); // A11 * B11

    return totalFlops;
}

MatrixXd decompressMatrix(const SVDNode* node) {
    if (!node) return MatrixXd();

    if (node->isLeaf) {
        // Rekonstrukcja bloku: A ~ U * S * V^T
        // Jeśli S jest puste (macierz zerowa), zwracamy zera o odpowiednich wymiarach
        if (node->S.size() == 0) {
            return MatrixXd::Zero(node->U.rows(), node->V.rows());
        }
        return node->U * node->S.asDiagonal() * node->V.transpose();
    } else {
        // Rekurencyjne pobranie 4 ćwiartek
        MatrixXd TL = decompressMatrix(node->children[0].get()); // Top-Left
        MatrixXd TR = decompressMatrix(node->children[1].get()); // Top-Right
        MatrixXd BL = decompressMatrix(node->children[2].get()); // Bottom-Left
        MatrixXd BR = decompressMatrix(node->children[3].get()); // Bottom-Right

        // Scalenie w jedną dużą macierz
        MatrixXd result(TL.rows() + BL.rows(), TL.cols() + TR.cols());
        
        // Operator przecinkowy Eigen do składania bloków
        result << TL, TR,
                  BL, BR;
                  
        return result;
    }
}