#pragma once

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Spectra/SymEigsSolver.h>
#include <memory>

using namespace Eigen;
using namespace std;
using namespace Spectra;

struct SVDNode {
    bool isLeaf;
    MatrixXd U;
    VectorXd S;
    MatrixXd V;
    std::vector<std::unique_ptr<SVDNode>> children;

    SVDNode() : isLeaf(true) {}
};

void computePartialSVD(const MatrixXd& A, int k, MatrixXd& U_out, VectorXd& S_out, MatrixXd& V_out);