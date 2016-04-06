#include <Eigen/Dense>
#include <iostream>

//! \brief Given a matrix A of linearly independent columns, returns Gram-Schmidt orthonormalization
//! Ustable GS algorithm. Output is prone to cancellation issues.
//! \param[in] A Matrix of linearly independent columns
//! \return Matrix with ONB of $span(a_1, \cdots, a_n)$
template <class Matrix>
Matrix gramschmidt( const Matrix & A ) {
    
    Matrix Q = A;
    // First vector just gets normalized
    Q.col(0).normalize();
    
    for(unsigned int j = 1; j < A.cols(); ++j) {
        // Replace inner loop over each previous vector in Q with fast matrix-vector multiplication
        Q.col(j) -= Q.leftCols(j) * (Q.leftCols(j).transpose() * A.col(j));
        
        // Normalize vector if possible (othw. means colums of A almsost lin. dep.
        if( Q.col(j).norm() <= 10e-14 * A.col(j).norm() ) {
            std::cerr << "Gram-Schmidt failed because A has lin. dep columns. Bye." << std::endl;
            break;
        } else {
            Q.col(j).normalize();
        }
    }
    
    return Q;
}

int main(void) {
    // Ortho test
    unsigned int n = 9;
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
    Eigen::MatrixXd Q = gramschmidt( A );
    
    // Output should be idenity matrix
    std::cout << Q*Q.transpose() << std::endl;
}
