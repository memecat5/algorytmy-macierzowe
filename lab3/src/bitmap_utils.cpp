
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

void saveSeparateChannels(const string& prefix, const MatrixXd& R, const MatrixXd& G, const MatrixXd& B) {
    int rows = R.rows();
    int cols = R.cols();
    
    // Macierz samych zer (czarny kolor)
    MatrixXd Zeros = MatrixXd::Zero(rows, cols);

    // 1. Zapisz tylko kanał R (R, 0, 0)
    savePPM(prefix + "_red.ppm", R, Zeros, Zeros);
    
    // 2. Zapisz tylko kanał G (0, G, 0)
    savePPM(prefix + "_green.ppm", Zeros, G, Zeros);

    // 3. Zapisz tylko kanał B (0, 0, B)
    savePPM(prefix + "_blue.ppm", Zeros, Zeros, B);
    
    // Opcjonalnie: Zapisz jako skala szarości (często łatwiej tak analizować detale)
    // savePPM(prefix + "_red_mono.ppm", R, R, R);
    // savePPM(prefix + "_green_mono.ppm", G, G, G);
    // savePPM(prefix + "_blue_mono.ppm", B, B, B);

    cout << "Zapisano osobne kanaly: " << prefix << "_red/green/blue.ppm" << endl;
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

void drawRankMap(const SVDNode* node, Eigen::Ref<Eigen::MatrixXd> target, int max_b) {
    if (node->isLeaf) {
        // 1. Obliczamy intensywność koloru (0 - 255)
        // Ranga 0/1 -> Ciemny szary
        // Ranga max_b -> Biały
        int currentRank = node->S.size();
        
        // Zabezpieczenie przez dzieleniem przez 0
        if (max_b < 1) max_b = 1;

        double intensity = (double)currentRank / (double)max_b * 255.0;
        
        // Przycinamy do zakresu 0-255 (dla bezpieczeństwa)
        intensity = std::clamp(intensity, 0.0, 255.0);

        // 2. Wypełniamy cały blok tym kolorem
        target.setConstant(intensity);

        // 3. Rysujemy czarną ramkę (żeby było widać podział na kwadraty)
        // Rysujemy tylko jeśli blok jest większy niż 2x2 (żeby nie zamazać wszystkiego ramką)
        // if (target.rows() > 2 && target.cols() > 2) {
        //     double borderColor = 0.0; // Czarny
            
        //     // Górna i dolna krawędź
        //     target.row(0).setConstant(borderColor);
        //     target.row(target.rows() - 1).setConstant(borderColor);
            
        //     // Lewa i prawa krawędź
        //     target.col(0).setConstant(borderColor);
        //     target.col(target.cols() - 1).setConstant(borderColor);
        // }

    } else {
        // Rekurencja (standardowy podział Quadtree)
        int r = target.rows();
        int c = target.cols();
        int hR = r / 2;
        int hC = c / 2;

        drawRankMap(node->children[0].get(), target.block(0, 0, hR, hC), max_b);
        drawRankMap(node->children[1].get(), target.block(0, hC, hR, c - hC), max_b);
        drawRankMap(node->children[2].get(), target.block(hR, 0, r - hR, hC), max_b);
        drawRankMap(node->children[3].get(), target.block(hR, hC, r - hR, c - hC), max_b);
    }

    if (target.rows() > 2 && target.cols() > 2) {
            double borderColor = 0.0; // Czarny
            
            // Górna i dolna krawędź
            target.row(0).setConstant(borderColor);
            target.row(target.rows() - 1).setConstant(borderColor);
            
            // Lewa i prawa krawędź
            target.col(0).setConstant(borderColor);
            target.col(target.cols() - 1).setConstant(borderColor);
        }
}


// ---------------------- HSV ------------------------------
void RGBtoHSV(const MatrixXd& R, const MatrixXd& G, const MatrixXd& B,
              MatrixXd& H, MatrixXd& S, MatrixXd& V) {
    int rows = R.rows();
    int cols = R.cols();
    H.resize(rows, cols);
    S.resize(rows, cols);
    V.resize(rows, cols);

    for (int y = 0; y < rows; ++y) {
        for (int x = 0; x < cols; ++x) {
            double r = R(y, x);
            double g = G(y, x);
            double b = B(y, x);

            double maxVal = std::max({r, g, b});
            double minVal = std::min({r, g, b});
            double delta = maxVal - minVal;

            // --- Value (0 - 255) ---
            V(y, x) = maxVal;

            // --- Saturation (0.0 - 1.0) ---
            if (maxVal > 0.0) {
                S(y, x) = delta / maxVal;
            } else {
                S(y, x) = 0.0;
            }

            // --- Hue (0.0 - 360.0) ---
            if (delta < 0.00001) {
                H(y, x) = 0.0; // Szary, Hue nieistotne
            } else {
                if (r >= maxVal) {
                    H(y, x) = (g - b) / delta;
                    if (H(y, x) < 0.0) H(y, x) += 6.0;
                } else if (g >= maxVal) {
                    H(y, x) = 2.0 + (b - r) / delta;
                } else {
                    H(y, x) = 4.0 + (r - g) / delta;
                }
                H(y, x) *= 60.0;
            }
        }
    }
}

void HSVtoRGB(const MatrixXd& H, const MatrixXd& S, const MatrixXd& V,
              MatrixXd& R, MatrixXd& G, MatrixXd& B) {
    int rows = H.rows();
    int cols = H.cols();
    R.resize(rows, cols);
    G.resize(rows, cols);
    B.resize(rows, cols);

    for (int y = 0; y < rows; ++y) {
        for (int x = 0; x < cols; ++x) {
            double h = H(y, x);
            double s = S(y, x);
            double v = V(y, x);

            if (s <= 0.0) {
                // Skala szarości
                R(y, x) = v;
                G(y, x) = v;
                B(y, x) = v;
            } else {
                if (h >= 360.0) h = 0.0;
                h /= 60.0;
                long i = (long)h;
                double ff = h - i;
                double p = v * (1.0 - s);
                double q = v * (1.0 - (s * ff));
                double t = v * (1.0 - (s * (1.0 - ff)));

                switch (i) {
                case 0: R(y, x) = v; G(y, x) = t; B(y, x) = p; break;
                case 1: R(y, x) = q; G(y, x) = v; B(y, x) = p; break;
                case 2: R(y, x) = p; G(y, x) = v; B(y, x) = t; break;
                case 3: R(y, x) = p; G(y, x) = q; B(y, x) = v; break;
                case 4: R(y, x) = t; G(y, x) = p; B(y, x) = v; break;
                default: R(y, x) = v; G(y, x) = p; B(y, x) = q; break;
                }
            }
        }
    }
}

void saveHSVChannels(const std::string& prefix, const MatrixXd& H, const MatrixXd& S, const MatrixXd& V) {
    // Żeby zapisać H i S jako obrazki PPM, musimy je przeskalować do 0-255
    
    // H: 0-360 -> 0-255
    MatrixXd H_vis = H * (255.0 / 360.0);
    
    // S: 0-1 -> 0-255
    MatrixXd S_vis = S * 255.0;
    
    // V jest już 0-255
    
    savePPM(prefix + "_channel_H.ppm", H_vis, H_vis, H_vis);
    savePPM(prefix + "_channel_S.ppm", S_vis, S_vis, S_vis);
    savePPM(prefix + "_channel_V.ppm", V, V, V);
    
    std::cout << "Zapisano wizualizacje kanalow HSV." << std::endl;
}