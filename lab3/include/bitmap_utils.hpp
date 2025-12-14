#include <Eigen/Core>

using namespace std;
using namespace Eigen;

void savePPM(const string& filename, const MatrixXd& R, const MatrixXd& G, const MatrixXd& B);
void loadPPM(const string& filename, MatrixXd& R, MatrixXd& G, MatrixXd& B);
void generateBitmap(int w, int h, MatrixXd& R, MatrixXd& G, MatrixXd& B);
