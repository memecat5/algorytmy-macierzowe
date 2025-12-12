#pragma once
#include <Eigen/Core>
#include <Eigen/SVD>
#include <memory>

#include "partial_svd.hpp"

std::unique_ptr<SVDNode> compressMatrixRecursive(const MatrixXd& block, double delta, int b);