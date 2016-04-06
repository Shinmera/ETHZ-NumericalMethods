#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <iostream>

//! \brief Compute y = A*x for A banded matrix with diagonal structure
//! \param[in] a n-1 dim. vector for upper diagonal
//! \param[in] b n-2 dim. vector for second lower diagonal
//! \param[in] x n dim. vector for Ax = y
//! \param[out] y n dim. vector y = Ax
template <class Vector>
void multAx(const Vector & a, const Vector & b, const Vector & x, Vector & y) {
    unsigned int n = x.size();
    if( a.size() < n - 1 || b.size() < n - 2 ) {
        // Size check, error if do not match
        std::cerr << "Error: size mismatch!" << std::endl;
        return;
    }
    y.resize(n); // Ensure y has size n
    
    // Handle first two rows
    if(n > 1) y(0) = 2*x(0) + a(0)*x(1);
    else { // Special case n = 1
        y(0) = 2*x(0);
        return;
    }
    if(n > 2) y(1) = 2*x(1) + a(1)*x(2);
    else { // Special case n = 2
        y(1) = 2*x(1);
        return;
    }
    
    // Row if n > 3 and without last row
    for(unsigned int i = 2; i < n-1; ++i) {
        y(i) = 2*x(i) + b(i-2)*x(i-2) + a(i)*x(i+1);
    }
    
    // Last row special case
    if(n > 2) y(n-1) = 2*x(n-1) + b(n-3)*x(n-3);
}

//! \brief Solve y = A*x for A banded matrix with upper tirangular sparse structure
//! \param[in] a n-1 dim. vector for upper diagonal
//! \param[in] r n dim. vector for Ax = r
//! \param[out] x n dim. solution vector r = Ax
template <class Vector>
void solvelseAupper(const Vector & a, const Vector & r, Vector & x) {
    // Dimensions + ensure correct dimensions
    unsigned int n = r.size();
    x.resize(n);
    
    // Backward substitution
    x(n-1) = 0.5 * r(n-1);
    for(int j = n-2; j >= 0; --j) {
        x(j) = 0.5*(r(j) - a(j)*x(j+1));
    }
}

//! \brief Solve y = A*x for A banded matrix using Gaussian-elimination (no pivot)
//! \param[in] a n-1 dim. vector for upper diagonal
//! \param[in] b n-2 dim. vector for second lower diagonal
//! \param[in] r n dim. vector for Ax = r
//! \param[out] x n dim. solution vector r = Ax
template <class Vector>
void solvelseA(const Vector & a, const Vector & b, const Vector & r, Vector & x) {
    // Setup type, size and copy data
    typedef typename Vector::Scalar Scalar;
    int n = r.size();
    x = r;
    
    // Fill in matrix, we reserve 5 nnz per row for Gauss fill-in
    Eigen::SparseMatrix<Scalar> A(n,n);
    A.reserve(5);
    for(int i = 0; i < n; ++i) {
        A.insert(i,i) = 2;
        if(i < n-1) A.insert(i,i+1) = a(i);
        if(i >= 2)  A.insert(i,i-2) = b(i-2);
    }
    A.makeCompressed();
    
    // 1st stage: elimination (in place)
    for(int i = 0; i < n-1; ++i) {
        for(int k = i+1; k < std::min(i+3,n); ++k) {
            Scalar fac = A.coeffRef(k,i)/A.coeffRef(i,i);
            for(int l = i; l < std::min(i+3,n); ++l) {
                A.coeffRef(k,l) -= fac * A.coeffRef(i,l);
            }
            x(k) -= fac * x(i);
        }
    }
    
    // 2nd stage: backwards substitution (in place)
    x(n-1) /= A.coeffRef(n-1,n-1);
    for(int i = n-2; i >= 0; --i) {
        for(int k = i+1; k < std::min(i+3,n); ++k) {
            x(i) -= x(k)*A.coeffRef(i,k);
        }
        x(i) /= A.coeffRef(i,i);
    }
}

//! \brief Solve y = A*x for A banded matrix using SparseLU
//! \param[in] a n-1 dim. vector for upper diagonal
//! \param[in] b n-2 dim. vector for second lower diagonal
//! \param[in] r n dim. vector for Ax = r
//! \param[out] x n dim. solution vector r = Ax
template <class Vector>
void solvelseAEigen(const Vector & a, const Vector & b, const Vector & r, Vector & x) {
    // Dimensions
    typedef typename Vector::Scalar Scalar;
    unsigned int n = r.size();
    
    // Fill in matrix + reserve
    Eigen::SparseMatrix<Scalar> A(n,n);
    A.reserve(3);
    for(unsigned int i = 0; i < n; ++i) {
        A.insert(i,i) = 2;
        if(i < n-1) A.insert(i,i+1) = a(i);
        if(i >= 2)  A.insert(i,i-2) = b(i-2);
    }
    A.makeCompressed();
    
    // Call SparseLU
    Eigen::SparseLU< Eigen::SparseMatrix<Scalar> >   solver;
    solver.analyzePattern(A); 
    solver.compute(A);
    x = solver.solve(r);
}

int main() {
    // Random vectors to check correctness
    std::cout << "*** Check correctness of  banded solvers implementations" << std::endl;
    int n = 9;
    Eigen::VectorXd a = Eigen::VectorXd::Random(n-1);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n-2);
    Eigen::VectorXd y = Eigen::VectorXd::Random(n);
    Eigen::VectorXd x;
    
    std::cout << "Original:" << y << std::endl;
    
    // Compute y = A*inv(A)*y
    solvelseAupper(a,y,x);
    multAx(a,b,x,y);
    
    b = Eigen::VectorXd::Random(n-2);
    
    // Should be same as before
    std::cout << "Upper:" << y << std::endl;
    
    // Compute y = A*inv(A)*y
    solvelseA(a,b,y,x);
    multAx(a,b,x,y);
    
    // Should be same as before
    std::cout << "Own:" << y << std::endl;
    
    // Compute y = A*inv(A)*y
    solvelseAEigen(a,b,y,x);
    multAx(a,b,x,y);
    
    // Should be same as before
    std::cout << "Eigen:" << y << std::endl;
    
}
