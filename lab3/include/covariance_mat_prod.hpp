#pragma once

#include <Spectra/SymEigsSolver.h>

using namespace Eigen;
using namespace std;
using namespace Spectra;

// ==========================================
// Klasa Operatora dla Spectra (A^T * A)
// ==========================================
class CovarianceMatProd {
private:
    const MatrixXd& mat_A;

public:
    // Spectra potrzebuje wiedzieć, jaki typ danych jest używany
    using Scalar = double;

    CovarianceMatProd(const MatrixXd& A) : mat_A(A) {}

    int rows() const { return mat_A.cols(); }
    int cols() const { return mat_A.cols(); }

    // y_out = (A^T * A) * x_in
    void perform_op(const double* x_in, double* y_out) const {
        Map<const VectorXd> x(x_in, mat_A.cols());
        Map<VectorXd> y(y_out, mat_A.cols());

        // Krok 1: temp = A * x
        VectorXd temp = mat_A * x;

        // Krok 2: y = A^T * temp
        y = mat_A.transpose() * temp;
    }
};