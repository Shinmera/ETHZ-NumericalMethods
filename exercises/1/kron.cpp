#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "timer.h"

//! \brief Compute the Kronecker product $C = A \otimes B$.
//! \param[in] A Matrix $n \times n$
//! \param[in] B Matrix $n \times n$
//! \param[out] C Kronecker product of A and B of dim $n^2 \times n^2$
template <class Matrix>
void kron(const Matrix & A, const Matrix & B, Matrix & C)
{
    // TODO
}

//! \brief Compute the Kronecker product C = $A \otimes B$. Exploit matrix-vector product.
//! A,B and x must have dimension n \times n resp. n
//! \param[in] A Matrix $n \times n$
//! \param[in] B Matrix $n \times n$
//! \param[in] x Vector of dim $n$
//! \param[out] y Vector y = kron(A,B)*x
template <class Matrix, class Vector>
void kron_fast(const Matrix & A, const Matrix & B, const Vector & x, Vector & y)
{
    // TODO
}

//! \brief Compute the Kronecker product $C = A \otimes B$. Uses fast remapping tecniques (similar to Matlab reshape)
//! A,B and x must have dimension n \times n resp. n*n
//! \param[in] A Matrix $n \times n$
//! \param[in] B Matrix $n \times n$
//! \param[in] x Vector of dim $n$
//! \param[out] y Vector y = kron(A,B)*x
template <class Matrix, class Vector>
void kron_super_fast(const Matrix & A, const Matrix & B, const Vector & x, Vector & y)
{
    // TODO
}


int main(void) {
    
    // Check if kron works, cf. 
    Eigen::MatrixXd A(2,2);
    A << 1, 2, 3, 4;
    Eigen::MatrixXd B(2,2);
    B << 5, 6, 7, 8;
    Eigen::MatrixXd C;
    
    Eigen::VectorXd x = Eigen::VectorXd::Random(4);
    Eigen::VectorXd y;
    kron(A,B,C);
    y = C*x;
    std::cout << "kron(A,B)=" << std::endl << C << std::endl;
    std::cout << "Using kron: y=       " << std::endl << y << std::endl;
    
    kron_fast(A,B,x,y);
    std::cout << "Using kron_fast: y=  " << std::endl << y << std::endl;
    kron_super_fast(A,B,x,y);
    std::cout << "Using kron_super_fast: y=  " << std::endl << y << std::endl;
    
    // Compute runtime of different implementations of kron
    unsigned int repeats = 10;
    timer<> tm_kron, tm_kron_fast, tm_kron_super_fast;
    std::vector<int> times_kron, times_kron_fast, times_kron_super_fast;
    
    for(unsigned int p = 2; p <= 10; p++) {
        tm_kron.reset();
        tm_kron_fast.reset();
        tm_kron_super_fast.reset();
        for(unsigned int r = 0; r < repeats; ++r) {
            unsigned int M = pow(2,p);
            A = Eigen::MatrixXd::Random(M,M);
            B = Eigen::MatrixXd::Random(M,M);
            x = Eigen::VectorXd::Random(M*M);
            
            // May be too slow for large p, comment if so
            tm_kron.start();
//             kron(A,B,C);
//             y = C*x;
            tm_kron.stop();
            
            tm_kron_fast.start();
            kron_fast(A,B,x,y);
            tm_kron_fast.stop();
            
            tm_kron_super_fast.start();
            kron_super_fast(A,B,x,y);
            tm_kron_super_fast.stop();
        }
        
        std::cout << "Lazy Kron took:       " << tm_kron.avg().count() / 1000000. << " ms" << std::endl;
        std::cout << "Kron fast took:       " << tm_kron_fast.avg().count() / 1000000. << " ms" << std::endl;
        std::cout << "Kron super fast took: " << tm_kron_super_fast.avg().count() / 1000000. << " ms" << std::endl;
        times_kron.push_back( tm_kron.avg().count() );
        times_kron_fast.push_back( tm_kron_fast.avg().count() );
        times_kron_super_fast.push_back( tm_kron_super_fast.avg().count() );
    }
    
    for(auto it = times_kron.begin(); it != times_kron.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;
    for(auto it = times_kron_fast.begin(); it != times_kron_fast.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;
    for(auto it = times_kron_super_fast.begin(); it != times_kron_super_fast.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;
}
