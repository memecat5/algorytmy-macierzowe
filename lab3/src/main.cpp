#include <iostream>
#include <vector>
#include <string>
#include <memory>

// Biblioteki Eigen i Spectra
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Spectra/SymEigsSolver.h>

#include "partial_svd.hpp"
#include "compress.hpp"
#include "bitmap_utils.hpp"

using namespace Eigen;
using namespace std;
using namespace Spectra;


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

int main() {
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