#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <math.h>

//! \brief Compute powers of a square matrix using smart binary representation
//! \param[in,out] A matrix for which you want to compute $A^k$. $A^k$ is stored in $A$
//! \param[out] k integer for $A^k$
template <class Matrix>
void matPow(Matrix & A, unsigned int k) {
    // TODO
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

