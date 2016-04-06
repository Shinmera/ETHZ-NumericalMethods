#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <math.h>

//! \brief Compute powers of a square matrix using smart binary representation
//! \param[in,out] A matrix for which you want to compute $A^k$. $A^k$ is stored in $A$
//! \param[out] k integer for $A^k$
template <class Matrix>
void matPow(Matrix & A, unsigned int k) {
    Matrix X = Matrix::Identity(A.rows(), A.cols());
    
    // p is used as binary mask to check wether $k = \sum_{i = 0}^M b_i 2^i$ has 1 in the $i$-th binary digit
    // obtaining the binay representation of p can be done in many ways, here we use ~k & p to check i-th binary is 1
    unsigned int p = 1;
    // Cycle all the way up to the length of the binary representation of $k$
    for(unsigned int j = 1; j <= ceil(log2(k)); ++j) {
        if( ( ~k & p ) == 0 ) {
            X = X*A;
        }
        
        A = A*A;
        p = p << 1;
    }
    A = X;
}

int main(void) {
    // Check/Test with provided, complex, matrix
    unsigned int n = 3; // size of matrix
    unsigned int k = 9; // power
    
    double PI = M_PI; // from math.h
    std::complex<double> I = std::complex<double>(0,1); // imaginary unit
    
    Eigen::MatrixXcd A(n,n);
    
    for(unsigned int i = 0; i < n; ++i) {
        for(unsigned int j = 0; j < n; ++j) {
            A(i,j) = exp(2. * PI * I * (double) i * (double) j / (double) n) / sqrt((double) n);
        }
    }
    
    // Test with simple matrix/simple power
//     Eigen::MatrixXd A(2,2);
//     k = 3;
//     A << 1,2,3,4;
         
    // Output results
    std::cout << "A = " << A << std::endl;
    std::cout << "Eigen:" << std::endl << "A^" << k << " = " << A.pow(k) << std::endl;
    matPow(A, k);
    std::cout << "Ours:" << std::endl << "A^" << k << " = " << A <<std::endl;
}

