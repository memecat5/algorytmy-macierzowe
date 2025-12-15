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

VectorXd getSingularValues(const MatrixXd& M){
    BDCSVD<MatrixXd> svd(M, ComputeThinU | ComputeThinV);
    return svd.singularValues();
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

    VectorXd singular_values_R = getSingularValues(R);
    VectorXd singular_values_G = getSingularValues(G);
    VectorXd singular_values_B = getSingularValues(B);

    int b;
    int nth_delta;

    cout << "Podaj parametry dokladnosci - delta (po ktorej wartosci wlasnej obcinamy): ";
    cin >> nth_delta;

    if (nth_delta < 1 || nth_delta > 512){
        cout << "delta musi być indeksem wartości własnej\n";
        return -1;
    }

    cout << "b: ";
    cin >> b;

    cout << "Start b = " << b << "\ndeltaR = " <<
    singular_values_R[nth_delta-1] << "\tdeltaG = " << singular_values_G[nth_delta-1]
    << "\tdeltaB = " << singular_values_B[nth_delta-1] << endl;



    vector<std::unique_ptr<SVDNode>> channels(3);
    channels[0] = compressMatrixRecursive(R, singular_values_R[nth_delta-1], b);
    channels[1] = compressMatrixRecursive(G, singular_values_G[nth_delta-1], b);
    channels[2] = compressMatrixRecursive(B, singular_values_B[nth_delta-1], b);

    MatrixXd oR(H,W), oG(H,W), oB(H,W);
    drawCompressedMatrix(channels[0].get(), oR);
    drawCompressedMatrix(channels[1].get(), oG);
    drawCompressedMatrix(channels[2].get(), oB);

    saveSeparateChannels("final", oR, oG, oB);

    savePPM("final.ppm", oR, oG, oB);
    cout << "Zapisano final.ppm" << endl;
    return 0;
}