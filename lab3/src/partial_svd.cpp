#include <Eigen/Core>
#include <Eigen/SVD>
#include <Spectra/SymEigsSolver.h>

#include "covariance_mat_prod.hpp"

using namespace Eigen;
using namespace std;
using namespace Spectra;


/*
    Funkcja licząca PartialSVD, oparta na funkcji ze Spectry,
    czyli nakładki na Eigen. Dla odpowiednio małych macierzy
    liczy zwykłe SVD
*/
void computePartialSVD(const MatrixXd& A, int k, MatrixXd& U_out, VectorXd& S_out, MatrixXd& V_out) {
    int m = A.rows();
    int n = A.cols();

    int ncv = std::min({m, n, 2 * k + 5}); // Bezpieczny margines

    // WARUNEK KRYTYCZNY:
    // 1. Macierz musi być dostatecznie duża (np. > 16x16), żeby opłacało się uruchamiać Arnoldiego.
    // 2. ncv musi być ściśle większe od k (wymóg matematyczny algorytmu).
    // 3. ncv musi być mniejsze lub równe wymiarowi macierzy.
    bool useSpectra = (m > 16) && (n > 16) && (ncv > k) && (ncv <= std::min(m, n));

    if (useSpectra) {
        // operator dla Spectry
        CovarianceMatProd op(A);

        // konfiguracja solvera
        SymEigsSolver<CovarianceMatProd> eigs(op, k, ncv);

        eigs.init();
        
        int nconv = eigs.compute(Spectra::SortRule::LargestMagn, 1000, 1e-10);

        if (eigs.info() == CompInfo::Successful) {
            VectorXd eigenValues = eigs.eigenvalues();
            MatrixXd eigenVectors = eigs.eigenvectors();

            S_out.resize(k);
            for(int i=0; i<k; ++i) {
                // Zabezpieczenie przed ujemnym zerem numerycznym
                S_out(i) = std::sqrt(std::max(0.0, eigenValues(i)));
            }

            V_out = eigenVectors;

            U_out.resize(m, k);
            for(int i=0; i<k; ++i) {
                if (S_out(i) > 1e-9) {
                    U_out.col(i) = (A * V_out.col(i)) / S_out(i);
                } else {
                    U_out.col(i).setZero();
                }
            }
            return;
        }
    }

    // zwykłe SVD dla małych macierzy lub kiedy obliczenia się wywalą
    BDCSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
    int actualK = std::min((int)svd.singularValues().size(), k);
    U_out = svd.matrixU().leftCols(actualK);
    S_out = svd.singularValues().head(actualK);
    V_out = svd.matrixV().leftCols(actualK);
}