#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <memory>
#include <fstream>

// Biblioteki Eigen i Spectra
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Spectra/SymEigsSolver.h>

#include "covariance_mat_prod.hpp"
#include "partial_svd.hpp"
#include "compress.hpp"

using namespace Eigen;
using namespace std;
using namespace Spectra;

// ==========================================
// Funkcje pomocnicze
// ==========================================

void savePPM(const string& filename, const MatrixXd& R, const MatrixXd& G, const MatrixXd& B) {
    ofstream file(filename, ios::binary);
    if (!file.is_open()) return;
    file << "P6\n" << R.cols() << " " << R.rows() << "\n255\n";
    for (int y = 0; y < R.rows(); ++y) {
        for (int x = 0; x < R.cols(); ++x) {
            file << (unsigned char)std::clamp(R(y,x),0.0,255.0)
                 << (unsigned char)std::clamp(G(y,x),0.0,255.0)
                 << (unsigned char)std::clamp(B(y,x),0.0,255.0);
        }
    }
}

void generateBitmap(int w, int h, MatrixXd& R, MatrixXd& G, MatrixXd& B) {
    R.resize(h, w); G.resize(h, w); B.resize(h, w);
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            double v = std::sin(x*0.05) * std::cos(y*0.05);
            R(y, x) = 127 + 127 * v;
            G(y, x) = 127 + 127 * std::sin((x+y)*0.03);
            B(y, x) = (x % 50 < 25 && y % 50 < 25) ? 200 : 50; 
        }
    }
}

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
    generateBitmap(W, H, R, G, B);
    savePPM("original.ppm", R, G, B);

    int b = 4;        
    double delta = 15.0; 

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