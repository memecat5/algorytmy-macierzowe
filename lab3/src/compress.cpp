#include "compress.hpp"
#include "partial_svd.hpp"


/*
    Funkcja odpowiedzialna za rekurenycjną kompresję macierzy.
    Korzysta z computePartialSVD, czyli liczy tylko częściowe svd
    (k największych wartości osobliwych), ale dla odpowiednio małych
    macierzy schodzi do zwykłego SVD z Eigen'a (bo będzie szybciej)
*/
std::unique_ptr<SVDNode> compressMatrixRecursive(const MatrixXd& block, double delta, int b) {
    auto node = std::make_unique<SVDNode>();
    int rows = block.rows();
    int cols = block.cols();

    
    // taka optymalizacja jest już w computePartialSVD, więc tu chyba niepotrzebna
    // jednak była potrzebna XDD
    if (rows <= 2 || cols <= 2) {
        BDCSVD<MatrixXd> svd(block, ComputeThinU | ComputeThinV);
        node->isLeaf = true;
        node->U = svd.matrixU();
        node->S = svd.singularValues();
        node->V = svd.matrixV();
        return node;
    }

    MatrixXd U_part, V_part;
    VectorXd S_part;
    computePartialSVD(block, b, U_part, S_part, V_part);

    // najmniejsza wartość osobliwa z SVD
    double smallestCalculatedSigma = (S_part.size() > 0) ? S_part(S_part.size() - 1) : 0.0;
    
    // Sprawdzenie czy przybliżenie jest wystarczające
    bool goodApproximation = (S_part.size() < b) || (smallestCalculatedSigma < delta);

    if (goodApproximation) {
        // Zapisanie skompresowanej macierzy jako liść
        node->isLeaf = true;
        int keep = 0;
        for(int i = 0; i < S_part.size(); ++i) {
            if (S_part(i) >= delta) keep++;
        }
        if (keep == 0) keep = 1;

        node->U = U_part.leftCols(keep);
        node->S = S_part.head(keep);
        node->V = V_part.leftCols(keep);
    } else {
        // Podział macierzy na 4 i zejście niżej
        node->isLeaf = false;
        int hRows = rows / 2;
        int hCols = cols / 2;
        node->children.push_back(compressMatrixRecursive(block.block(0, 0, hRows, hCols), delta, b));
        node->children.push_back(compressMatrixRecursive(block.block(0, hCols, hRows, cols - hCols), delta, b));
        node->children.push_back(compressMatrixRecursive(block.block(hRows, 0, rows - hRows, hCols), delta, b));
        node->children.push_back(compressMatrixRecursive(block.block(hRows, hCols, rows - hRows, cols - hCols), delta, b));
    }

    return node;
}