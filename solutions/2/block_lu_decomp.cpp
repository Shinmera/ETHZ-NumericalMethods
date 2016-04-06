#include <Eigen/Dense>
#include <iostream>

#include "timer.h"

//! \brief Use efficient implementation A*x = bb
//! \param[in] R Matrix is nxn and upper triangular
//! \param[in] v Vector is nx1
//! \param[in] u Vector is nx1
//! \param[in] bb vector is (n+1)x1 and is stacked (b, \beta)^T =: b
//! \param[out] x solution A*bb = x
template <class Matrix, class Vector>
void solvelse(const Matrix & R, const Vector & v, const Vector & u,
              const Vector & bb, Vector & x) {
    // Size of R, which is wize of u, v, and size of bb is n+1
    unsigned int n = R.rows();
    if( n != R.cols() || n != u.size() || n != v.size() || n+1  != bb.size() ) {
        // Size check, error if do not match
        std::cerr << "Error: size mismatch!" << std::endl;
        return;
    }
    
    // This gets the scalar type of the matrix
    typedef typename Matrix::Scalar Scalar; 
    // Tell eigen R is triangular, we use the construct R.template to avoid compilation errors
    // cf. http://eigen.tuxfamily.org/dox-devel/TopicTemplateKeyword.html
    // here auto = const Eigen::TriangularView<const Matrix, Eigen::Upper>
    auto triR = R.template triangularView<Eigen::Upper>();
    
    // s is the Schur's complement and, in this case, is a scalar (and so is b_s)
    // snv = s^{-1}, b_s as in lecture notes
    // sinvbs = s^{-1}*b_s
    Scalar sinv = - 1. / u.dot(triR.solve(v));
    Scalar bs = (bb(n) - u.dot(triR.solve(bb.head(n))));
    Scalar sinvbs = sinv*bs;
    
    // Stack the vector (z, \xi)^T =: x
    x = Vector::Zero(n+1);
    x << triR.solve(bb.head(n) - v*sinvbs), sinvbs;
}

//! \brief Use Eigen's LU-solver to solve Ax = y
//! \param[in] R Matrix is nxn and upper triangular
//! \param[in] v Vector is nx1
//! \param[in] u Vector is nx1
//! \param[in] bb vector is (n+1)x1 and is stacked (b, \beta)^T =: b
//! \param[out] x solution A*bb = x
template <class Matrix, class Vector>
void solvelse_lu(const Matrix & R, const Vector & v, const Vector & u,
                 const Vector & bb, Vector & x) {
    // Size of R, which is wize of u, v, and size of bb is n+1
    unsigned int n = R.rows();
    if( n != R.cols() && n != u.size() && n != v.size() && n+1  != bb.size() ) {
        // Size check, error if do not match
        std::cerr << "Error: size mismatch!" << std::endl;
        return;
    }
    
    Matrix A(n+1,n+1);
    A << R, v, u.transpose(), 0;
    
    x = A.partialPivLu().solve(bb);
}

int main() {
    // Bunch of random vectors/matrices
    int n = 9;
    Eigen::MatrixXd R = Eigen::MatrixXd::Random(n,n).triangularView<Eigen::Upper>();
    Eigen::VectorXd v = Eigen::VectorXd::Random(n);
    Eigen::VectorXd u = Eigen::VectorXd::Random(n);
    Eigen::VectorXd bb = Eigen::VectorXd::Random(n+1);
    Eigen::VectorXd x;
    
    // Check that answer is correct
    std::cout << "*** Check correctness of Block Gauss implementation" << std::endl;
    solvelse(R, v, u, bb, x);
    std::cout << "Block Gauss:" << std::endl << x << std::endl;
    // Compare with eigen partial pivot LU-solve
    solvelse_lu(R, v, u, bb, x);
    std::cout << "Eigen LU:" << std::endl << x << std::endl;
    
    // Compute runtime of different implementations of kron
    std::cout << "*** Runtime comparisons of Block Gauss implementation VS Eigen LU" << std::endl;
    unsigned int repeats = 3;
    timer<> tm_own, tm_eigen_lu;
    
    for(unsigned int p = 2; p <= 10; p++) {
        tm_own.reset();
        tm_eigen_lu.reset();
        for(unsigned int r = 0; r < repeats; ++r) {
            unsigned int M = pow(2,p);
            R = Eigen::MatrixXd::Random(M,M).triangularView<Eigen::Upper>(); // We will use only upper triangular part
            v = Eigen::VectorXd::Random(M);
            u = Eigen::VectorXd::Random(M);
            bb = Eigen::VectorXd::Random(M+1);
            
            tm_own.start();
            solvelse(R, v, u, bb, x);
            tm_own.stop();
            
            tm_eigen_lu.start();
            solvelse_lu(R, v, u, bb, x);
            tm_eigen_lu.stop();
        }
        
        std::cout << "Own implementation took:       " << tm_own.min().count() / 1000000. << " ms" << std::endl;
        std::cout << "Eigen solver took:             " << tm_eigen_lu.min().count() / 1000000. << " ms" << std::endl;
    }
}
