#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <fstream> // Do zapisu CSV
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/SVD>

using namespace Eigen;
using namespace std;

// --- STRUKTURY I ALGORYTMY H-MACIERZY (Z poprzednich zadań) ---

struct HNode {
    bool is_leaf;
    int t_min, t_max, s_min, s_max;
    int rank;
    int rows, cols;
    MatrixXd U, V;
    VectorXd SingularValues;
    std::vector<HNode*> children;

    HNode(int tmin, int tmax, int smin, int smax) 
        : t_min(tmin), t_max(tmax), s_min(smin), s_max(smax), is_leaf(false), rank(0) {
        rows = tmax - tmin + 1;
        cols = smax - smin + 1;
    }
    ~HNode() { for (auto c : children) delete c; }
};

MatrixXd GetBlock(const MatrixXd& A, int t1, int t2, int s1, int s2) {
    return A.block(t1, s1, t2 - t1 + 1, s2 - s1 + 1);
}

// Rekurencyjna kompresja [cite: 1191]
HNode* CompressMatrix(const MatrixXd& A, int t_min, int t_max, int s_min, int s_max, int max_rank, double epsilon) {
    HNode* node = new HNode(t_min, t_max, s_min, s_max);
    int rows = node->rows;
    int cols = node->cols;

    if (rows <= 32 || cols <= 32) { // Warunek stopu dla małych bloków
        node->is_leaf = true;
        BDCSVD<MatrixXd> svd(GetBlock(A, t_min, t_max, s_min, s_max), ComputeThinU | ComputeThinV);
        node->rank = svd.rank();
        node->U = svd.matrixU();
        node->SingularValues = svd.singularValues();
        node->V = svd.matrixV();
        return node;
    }

    MatrixXd block = GetBlock(A, t_min, t_max, s_min, s_max);
    BDCSVD<MatrixXd> svd(block, ComputeThinU | ComputeThinV);
    const VectorXd& vals = svd.singularValues();
    
    int k = 0;
    for(int i=0; i<vals.size(); ++i) if(vals(i) > epsilon) k++;

    // Warunek dopuszczalności
    if (k <= max_rank && k <= std::min(rows, cols) / 2) {
        node->is_leaf = true;
        node->rank = k;
        if (k > 0) {
            node->U = svd.matrixU().leftCols(k);
            node->SingularValues = vals.head(k);
            node->V = svd.matrixV().leftCols(k);
        } else {
             node->U = MatrixXd(rows, 0); node->V = MatrixXd(cols, 0);
        }
    } else {
        int t_mid = t_min + rows/2 - 1;
        int s_mid = s_min + cols/2 - 1;
        node->children.push_back(CompressMatrix(A, t_min, t_mid, s_min, s_mid, max_rank, epsilon));
        node->children.push_back(CompressMatrix(A, t_min, t_mid, s_mid+1, s_max, max_rank, epsilon));
        node->children.push_back(CompressMatrix(A, t_mid+1, t_max, s_min, s_mid, max_rank, epsilon));
        node->children.push_back(CompressMatrix(A, t_mid+1, t_max, s_mid+1, s_max, max_rank, epsilon));
    }
    return node;
}

// Dekonstrukcja do macierzy gęstej (do liczenia błędu) [cite: 1205, 1220]
MatrixXd Decompress(HNode* node) {
    if (node->is_leaf) {
        if (node->rank > 0) return node->U * node->SingularValues.asDiagonal() * node->V.transpose();
        return MatrixXd::Zero(node->rows, node->cols);
    }
    MatrixXd R(node->rows, node->cols);
    int h = node->children[0]->rows;
    int w = node->children[0]->cols;
    R.block(0, 0, h, w) = Decompress(node->children[0]);
    R.block(0, w, h, node->cols-w) = Decompress(node->children[1]);
    R.block(h, 0, node->rows-h, w) = Decompress(node->children[2]);
    R.block(h, w, node->rows-h, node->cols-w) = Decompress(node->children[3]);
    return R;
}

// Mnożenie H-Macierz x Wektor  [cite: 1193, 2219-2249]
VectorXd HMatrixVectorMult(HNode* node, const VectorXd& x) {
    if (node->is_leaf) {
        if (node->rank == 0) return VectorXd::Zero(node->rows);
        return node->U * (node->SingularValues.asDiagonal() * (node->V.transpose() * x));
    }
    int split = node->children[0]->cols;
    VectorXd x1 = x.head(split);
    VectorXd x2 = x.tail(x.size() - split);
    VectorXd y1 = HMatrixVectorMult(node->children[0], x1) + HMatrixVectorMult(node->children[1], x2);
    VectorXd y2 = HMatrixVectorMult(node->children[2], x1) + HMatrixVectorMult(node->children[3], x2);
    VectorXd y(y1.size() + y2.size());
    y << y1, y2;
    return y;
}

// Pomocnicze dodawanie H-Macierzy (uproszczone przez dekompresję liści)
HNode* HMatrixAdd(HNode* A, HNode* B, double eps, int mr) {
    HNode* res = new HNode(A->t_min, A->t_max, A->s_min, A->s_max);

    // ZMIANA: Używamy || (lub). Jeśli chociaż jeden jest liściem, 
    // nie możemy schodzić niżej, musimy wykonać operację na poziomie danych.
    if (A->is_leaf || B->is_leaf) {
        res->is_leaf = true;
        
        // Decompress obsługuje rekurencyjnie węzły, więc zadziała 
        // nawet jeśli np. A jest Węzłem (zdekompresuje całe poddrzewo do macierzy)
        MatrixXd dA = Decompress(A);
        MatrixXd dB = Decompress(B);
        
        MatrixXd sum = dA + dB;

        // Ponowna kompresja wyniku
        BDCSVD<MatrixXd> svd(sum, ComputeThinU | ComputeThinV);
        int k = 0; 
        for(int i=0; i<svd.singularValues().size(); ++i) 
            if(svd.singularValues()(i) > eps) k++;
        
        if(k > mr) k = mr;
        res->rank = k;

        if(k > 0) {
            res->U = svd.matrixU().leftCols(k);
            res->SingularValues = svd.singularValues().head(k);
            res->V = svd.matrixV().leftCols(k);
        } else {
            res->U = MatrixXd(A->rows, 0); 
            res->V = MatrixXd(A->cols, 0);
        }
    } else {
        // Oba są węzłami wewnętrznymi - możemy bezpiecznie rekurencyjnie wywołać dla dzieci
        for(int i=0; i<4; ++i) {
            res->children.push_back(HMatrixAdd(A->children[i], B->children[i], eps, mr));
        }
    }
    return res;
}

// Mnożenie H-Macierz x H-Macierz  [cite: 1194, 2365-2368]
HNode* HMatrixMatrixMult(HNode* A, HNode* B, double eps, int mr) {
    HNode* res = new HNode(A->t_min, A->t_max, B->s_min, B->s_max);

    // ZMIANA: Używamy || (lub). Jeśli którykolwiek jest liściem, 
    // wykonujemy mnożenie gęste (dekompresja -> mnożenie -> kompresja)
    if (A->is_leaf || B->is_leaf) {
        res->is_leaf = true;
        
        MatrixXd dA = Decompress(A);
        MatrixXd dB = Decompress(B);
        
        MatrixXd prod = dA * dB;

        BDCSVD<MatrixXd> svd(prod, ComputeThinU | ComputeThinV);
        int k = 0; 
        for(int i=0; i<svd.singularValues().size(); ++i) 
            if(svd.singularValues()(i) > eps) k++;
        
        if(k > mr) k = mr;
        res->rank = k;

        if(k > 0) {
            res->U = svd.matrixU().leftCols(k);
            res->SingularValues = svd.singularValues().head(k);
            res->V = svd.matrixV().leftCols(k);
        } else { 
            res->U = MatrixXd(res->rows, 0); 
            res->V = MatrixXd(res->cols, 0); 
        }
        return res;
    }

    // Oba są węzłami wewnętrznymi - standardowa rekurencja blokowa
    // C11 = A11*B11 + A12*B21
    HNode* t1 = HMatrixMatrixMult(A->children[0], B->children[0], eps, mr);
    HNode* t2 = HMatrixMatrixMult(A->children[1], B->children[2], eps, mr);
    res->children.push_back(HMatrixAdd(t1, t2, eps, mr)); 
    delete t1; delete t2;

    // C12 = A11*B12 + A12*B22
    HNode* t3 = HMatrixMatrixMult(A->children[0], B->children[1], eps, mr);
    HNode* t4 = HMatrixMatrixMult(A->children[1], B->children[3], eps, mr);
    res->children.push_back(HMatrixAdd(t3, t4, eps, mr)); 
    delete t3; delete t4;

    // C21 = A21*B11 + A22*B21
    HNode* t5 = HMatrixMatrixMult(A->children[2], B->children[0], eps, mr);
    HNode* t6 = HMatrixMatrixMult(A->children[3], B->children[2], eps, mr);
    res->children.push_back(HMatrixAdd(t5, t6, eps, mr)); 
    delete t5; delete t6;

    // C22 = A21*B12 + A22*B22
    HNode* t7 = HMatrixMatrixMult(A->children[2], B->children[1], eps, mr);
    HNode* t8 = HMatrixMatrixMult(A->children[3], B->children[3], eps, mr);
    res->children.push_back(HMatrixAdd(t7, t8, eps, mr)); 
    delete t7; delete t8;
    
    return res;
}

// --- GENERATOR 3D (Zadanie 4) [cite: 1189] ---
MatrixXd Generate3DTopology(int k) {
    int N_side = std::pow(2, k);
    int N = N_side * N_side * N_side;
    MatrixXd A = MatrixXd::Zero(N, N);
    
    auto idx = [&](int x, int y, int z) { return z * N_side * N_side + y * N_side + x; };
    
    std::mt19937 gen(1234); // Stałe ziarno dla powtarzalności
    std::uniform_real_distribution<> dis(0.1, 1.0);

    for(int z=0; z<N_side; ++z) {
        for(int y=0; y<N_side; ++y) {
            for(int x=0; x<N_side; ++x) {
                int r = idx(x,y,z);
                A(r, r) = dis(gen) + 5.0; // Przekątna
                
                int nbs[6][3] = {{x-1,y,z}, {x+1,y,z}, {x,y-1,z}, {x,y+1,z}, {x,y,z-1}, {x,y,z+1}};
                for(auto& nb : nbs) {
                    int nx=nb[0], ny=nb[1], nz=nb[2];
                    if(nx>=0 && nx<N_side && ny>=0 && ny<N_side && nz>=0 && nz<N_side) {
                        A(r, idx(nx,ny,nz)) = dis(gen);
                    }
                }
            }
        }
    }
    return A;
}

int main() {
    ofstream csv("results_task4.csv");
    csv << "k,N,Time_Vec_ms,Time_Mat_ms,Error_Vec,Error_Mat\n";
    
    int max_rank = 10;
    double eps_rel = 1e-3;

    cout << "=== RAPORT ZADANIE 4 ===" << endl;

    for (int k : {2, 3, 4}) { // k=2, 3, 4 [cite: 1202]
        cout << "\nPrzetwarzanie k=" << k << "..." << endl;
        
        // 1. Generowanie i Kompresja
        MatrixXd A = Generate3DTopology(k);
        int N = A.rows();
        double eps = A.norm() * eps_rel;
        
        HNode* hA = CompressMatrix(A, 0, N-1, 0, N-1, max_rank, eps);
        
        // 2. Mnożenie Macierz-Wektor [cite: 1202]
        VectorXd x = VectorXd::Random(N);
        auto t1 = chrono::high_resolution_clock::now();
        VectorXd yH = HMatrixVectorMult(hA, x);
        auto t2 = chrono::high_resolution_clock::now();
        double time_vec = chrono::duration_cast<chrono::microseconds>(t2 - t1).count() / 1000.0;

        // Obliczanie błędu Wektor [cite: 1206]
        VectorXd yRef = A * x;
        double err_vec = (yRef - yH).squaredNorm();

        // 3. Mnożenie Macierz-Macierz (A*A) [cite: 1217]
        auto t3 = chrono::high_resolution_clock::now();
        HNode* hA2 = HMatrixMatrixMult(hA, hA, eps, max_rank);
        auto t4 = chrono::high_resolution_clock::now();
        double time_mat = chrono::duration_cast<chrono::milliseconds>(t4 - t3).count();

        // Obliczanie błędu Macierz (Pomiń dla k=4 ze względu na czas A*A gęstego)
        double err_mat = 0.0;
        if (k < 4) {
            MatrixXd A2_ref = A * A;
            MatrixXd A2_H = Decompress(hA2);
            err_mat = (A2_ref - A2_H).squaredNorm(); // Frobenius squared [cite: 1221]
        } else {
            cout << "(Pominiecie obliczania bledu A^2 dla k=4 - zbyt dlugi czas dla macierzy gestej)" << endl;
            err_mat = -1.0; 
        }

        cout << "N=" << N << " | T_Vec=" << time_vec << "ms | T_Mat=" << time_mat << "ms" << endl;
        cout << "Err_Vec=" << err_vec << " | Err_Mat=" << err_mat << endl;

        // Zapis do CSV
        csv << k << "," << N << "," << time_vec << "," << time_mat << "," << err_vec << "," << err_mat << "\n";

        delete hA; delete hA2;
    }
    csv.close();
    cout << "\nDane zapisane do results_task4.csv. Uruchom skrypt Python." << endl;
    return 0;
}