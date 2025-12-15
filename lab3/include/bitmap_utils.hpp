#include <Eigen/Core>

#include "partial_svd.hpp"

using namespace std;
using namespace Eigen;

void savePPM(const string& filename, const MatrixXd& R, const MatrixXd& G, const MatrixXd& B);
void saveSeparateChannels(const string& prefix, const MatrixXd& R, const MatrixXd& G, const MatrixXd& B);
void loadPPM(const string& filename, MatrixXd& R, MatrixXd& G, MatrixXd& B);
void generateBitmap(int w, int h, MatrixXd& R, MatrixXd& G, MatrixXd& B);
void drawRankMap(const SVDNode* node, Eigen::Ref<Eigen::MatrixXd> target, int max_b);

// do hsv
// Konwersja przestrzeni barw
void RGBtoHSV(const Eigen::MatrixXd& R, const Eigen::MatrixXd& G, const MatrixXd& B,
              Eigen::MatrixXd& H, Eigen::MatrixXd& S, Eigen::MatrixXd& V);

void HSVtoRGB(const MatrixXd& H, const MatrixXd& S, const MatrixXd& V,
              MatrixXd& R, MatrixXd& G, MatrixXd& B);

// Zapisuje kanały HSV jako osobne obrazki (skalując H i S do 0-255 dla widoczności)
void saveHSVChannels(const std::string& prefix, const Eigen::MatrixXd& H, const Eigen::MatrixXd& S, const Eigen::MatrixXd& V);
