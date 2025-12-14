#ifndef COMPRESS_HPP
#define COMPRESS_HPP

#include "partial_svd.hpp"

#include <Eigen/Core>
#include <vector>
#include <memory>

using namespace Eigen;

// Funkcja kompresująca
std::unique_ptr<SVDNode> compressMatrixRecursive(const MatrixXd& block, double delta, int b);

// --- Funkcje analityczne ---

// Zwraca liczbę wierszy/kolumn reprezentowanych przez węzeł
int getRows(const SVDNode* node);
int getCols(const SVDNode* node);

// Oblicza liczbę operacji FLOPs dla mnożenia macierz-wektor (MVM)
long long getMVMFlops(const SVDNode* node);

// Oblicza liczbę operacji FLOPs dla mnożenia macierz-macierz (MMM): A * B
long long getMMMFlops(const SVDNode* nodeA, const SVDNode* nodeB);

#endif