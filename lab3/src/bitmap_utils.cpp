
#include <Eigen/Core>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "bitmap_utils.hpp"

using namespace std;
using namespace Eigen;

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

void loadPPM(const string& filename, MatrixXd& R, MatrixXd& G, MatrixXd& B) {
    ifstream file(filename, ios::binary);
    if (!file.is_open()) {
        cerr << "Nie mozna otworzyc pliku: " << filename << endl;
        exit(1);
    }

    string header;
    int w, h, maxVal;
    file >> header;
    if (header != "P6") {
        cerr << "To nie jest format P6!" << endl;
        exit(1);
    }
    file >> w >> h >> maxVal;
    file.ignore(256, '\n'); // Pomin jeden znak bialy po naglowku

    R.resize(h, w);
    G.resize(h, w);
    B.resize(h, w);

    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            unsigned char pixel[3];
            file.read(reinterpret_cast<char*>(pixel), 3);
            R(y, x) = pixel[0];
            G(y, x) = pixel[1];
            B(y, x) = pixel[2];
        }
    }
    cout << "Wczytano " << filename << " (" << w << "x" << h << ")" << endl;
}

/*
    Generowanie jakiejś bitmapy wymyślone przez czata
*/
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