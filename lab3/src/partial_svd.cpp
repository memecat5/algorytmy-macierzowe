#include <Eigen/Core>
#include <Eigen/SVD>
#include <Spectra/SymEigsSolver.h>

#include "covariance_mat_prod.hpp"

using namespace Eigen;
using namespace std;
using namespace Spectra;

// ==========================================
// Funkcja SVD oparta na Spectra::SymEigsSolver
// ==========================================

void computePartialSVD(const MatrixXd& A, int k, MatrixXd& U_out, VectorXd& S_out, MatrixXd& V_out) {
    int m = A.rows();
    int n = A.cols();

    // Spectra wymaga, by ncv (liczba wektorów Arnoldiego) spełniała k < ncv <= n
    int ncv = std::min(n, 2 * k + 5);

    // Warunek użycia Spectra (duże macierze)
    bool useSpectra = (n > 20) && (m > 20) && (ncv < n);

    if (useSpectra) {
        // 1. Definiujemy operator
        CovarianceMatProd op(A);

        // 2. Konfigurujemy solver.
        // POPRAWKA: W nowym Spectra podajemy TYLKO typ operatora w szablonie.
        SymEigsSolver<CovarianceMatProd> eigs(op, k, ncv);

        eigs.init();
        
        // POPRAWKA: Regułę sortowania (LargestMagn) podajemy tutaj, jako argument funkcji.
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

    // --- FALLBACK DO EIGEN ---
    BDCSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
    int actualK = std::min((int)svd.singularValues().size(), k);
    U_out = svd.matrixU().leftCols(actualK);
    S_out = svd.singularValues().head(actualK);
    V_out = svd.matrixV().leftCols(actualK);
}