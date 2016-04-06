#include <Eigen/Dense>
#include <Eigen/LU>

#include <iostream>
#include "timer.h"

//! \brief Circular shift (downwards) of b
//! \param[in,out] b Vector is nx1, shifted downwards
template <class Vector>
void shift(Vector & b) {
    typedef typename Vector::Scalar Scalar;
    int n = b.size();
    
    Scalar temp = b(n-1);
    for(int k = n-2; k >= 0; --k) {
        b(k+1) = b(k);
    }
    b(0) = temp;
}

//! \brief Compute X = inv(A)*[b_1,...,b_n], b_i = i-th cyclic shift of b
//! \param[in] A Matrix is nx
//! \param[in] b Vector is nx1
//! \param[out] X solution nxn matrix X = inv(A)*[b_1,...,b_n]
template <class Matrix, class Vector>
void solvpermb(const Matrix & A, Vector & b, Matrix & X) {
    // Size of b, which is size of A
    int n = b.size();
    if( n != A.cols() or n != A.rows() ) {
        // Size check, error if do not match
        std::cerr << "Error: size mismatch!" << std::endl;
        return;
    }
    X.resize(n,n);
    
    for(int l = 0; l < n; ++l) {
        X.col(l) = A.fullPivLu().solve(b);
        
        shift(b);
    }
}

//! \brief Compute X = inv(A)*[b_1,...,b_n], b_i = i-th cyclic shift of b in O(n^3)
//! \param[in] A Matrix is nx
//! \param[in] b Vector is nx1
//! \param[out] X solution nxn matrix X = inv(A)*[b_1,...,b_n]
template <class Matrix, class Vector>
void solvpermb_on3(const Matrix & A, Vector & b, Matrix & X) {
    // Size of b, which is size of A
    int n = b.size();
    if( n != A.cols() or n != A.rows() ) {
        // Size check, error if do not match
        std::cerr << "Error: size mismatch!" << std::endl;
        return;
    }
    X.resize(n,n);
    
    // Notice here we reuse the LU factorization
    Eigen::FullPivLU<Matrix> LU = A.fullPivLu();
    
    for(int l = 0; l < n; ++l) {
        X.col(l) = LU.solve(b);
        
        shift(b);
    }
}
    
int main() {
    unsigned int n = 9;
    // Compute with both LU and reuse LU
    std::cout << "*** Check correctness of permutation solver" << std::endl;
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
    Eigen::VectorXd b = Eigen::VectorXd::Random(n);
    Eigen::MatrixXd X;
    std::cout << "b = " << std::endl << b << std::endl;
    solvpermb(A,b,X);
    std::cout << "Direct porting from MALTAB: " << std::endl << X << std::endl;
    std::cout << "A*X = " << std::endl << A*X << std::endl;
    solvpermb_on3(A,b,X);
    std::cout << "Reusing LU: " << std::endl << X << std::endl;
    std::cout << "A*X = " << std::endl << A*X << std::endl;
    
    // Compute runtime of different implementations of kron
    std::cout << "*** Runtime comparisons of LU vs reuse LU" << std::endl;
    unsigned int repeats = 3;
    timer<> tm_lu, tm_reuse_lu;
    
    for(unsigned int p = 2; p <= 7; p++) {
        tm_lu.reset();
        tm_reuse_lu.reset();
        for(unsigned int r = 0; r < repeats; ++r) {
            unsigned int M = pow(2,p);
            A = Eigen::MatrixXd::Random(M,M);
            b = Eigen::VectorXd::Random(M);
            
            tm_lu.start();
            solvpermb(A,b,X);
            tm_lu.stop();
            
            tm_reuse_lu.start();
            solvpermb_on3(A,b,X);
            tm_reuse_lu.stop();
        }
        
        std::cout << "LU took:       " << tm_lu.avg().count() / 1000000. << " ms" << std::endl;
        std::cout << "Reuse LU took: " << tm_reuse_lu.avg().count() / 1000000. << " ms" << std::endl;
    }
}
