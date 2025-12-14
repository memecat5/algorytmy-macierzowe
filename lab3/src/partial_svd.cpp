#include <Eigen/Core>
#include <Eigen/SVD>
#include <Spectra/SymEigsSolver.h>
#include <iostream> // For std::cerr (optional logging)

#include "covariance_mat_prod.hpp"

using namespace Eigen;
using namespace std;
using namespace Spectra;

void computePartialSVD(const MatrixXd& A, int k, MatrixXd& U_out, VectorXd& S_out, MatrixXd& V_out) {
    int m = A.rows();
    int n = A.cols();

    // Spectra requires ncv to satisfy: k < ncv <= n
    // A standard heuristic is 2*k + 5 usually works well
    int ncv = std::min(n, 2 * k + 5);

    // FIX 1: Increase the threshold. 
    // Standard SVD (BDCSVD) is extremely fast for n < 100-200. 
    // Using Spectra for small/medium matrices (20-100) causes the stability issues you saw.
    bool useSpectra = (n > 128) && (m > 128) && (ncv < n) && (k < n);

    if (useSpectra) {
        // FIX 2: Wrap in try-catch to handle "TridiagEigen" failures gracefully.
        try {
            CovarianceMatProd op(A);
            SymEigsSolver<CovarianceMatProd> eigs(op, k, ncv);

            eigs.init();
            
            // The error was thrown here. If it fails, we catch it and fall back to BDCSVD.
            int nconv = eigs.compute(Spectra::SortRule::LargestMagn, 1000, 1e-10);

            if (eigs.info() == CompInfo::Successful) {
                VectorXd eigenValues = eigs.eigenvalues();
                MatrixXd eigenVectors = eigs.eigenvectors();

                S_out.resize(k);
                for(int i=0; i<k; ++i) {
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
                return; // Success, return early
            }
        } catch (...) {
            // If Spectra throws (e.g. TridiagEigen failure), ignore and proceed to fallback
        }
    }

    // --- FALLBACK: Standard Eigen SVD ---
    // Used for small matrices OR if Spectra fails
    BDCSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
    int actualK = std::min((int)svd.singularValues().size(), k);
    
    U_out = svd.matrixU().leftCols(actualK);
    S_out = svd.singularValues().head(actualK);
    V_out = svd.matrixV().leftCols(actualK);
}