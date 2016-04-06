#include <Eigen/Dense>
#include <iostream>
#include <functional> // for std::function
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

// New version 1.1: using sparse matrices, dropping general nonlinear solver, using error instead of residual

//! \brief Implements a single step of the fixed point iteration $x^{(k+1)} = A(x^{(k)})^{-1} * b$
//! \tparam func type of the lambda function implementing A(x)
//! \tparam Vector type for the vector b, x, x_new $\in \mathbf{R}^2$
//! \param[in] A lambda function implementing A(x)
//! \param[in] b rhs vector $b \in \mathbf{R}^n$
//! \param[in] x previous step $x^{(k)}$
//! \param[out] x_new next step $x^{(k+1)}$
template <class func, class Vector>
void fixed_point_step(func&& A, const Vector & b, const Vector & x, Vector & x_new) {
    // Next step
    auto T = A(x);
    Eigen::SparseLU<Eigen::SparseMatrix<double>> Ax_lu;
    Ax_lu.analyzePattern(T); 
    Ax_lu.factorize(T);
    x_new = Ax_lu.solve(b);
}

//! \brief Implements a single step of the Netwon iteration for $x^{(k+1)}$
//! Exploits Sherman-Morrison-Woodbury formula for fast inversion of rank-one modification of a matrix.
//! \tparam func type of the lambda function implementing A(x)
//! \tparam Vector type for the vector b, x, x_new $\in \mathbf{R}^2$
//! \param[in] A lambda function implementing A(x)
//! \param[in] b rhs vector $b \in \mathbf{R}^n$
//! \param[in] x previous step $x^{(k)}$
//! \param[out] x_new next step in Newton iteration $x^{(k+1)}$
template <class func, class Vector>
void newton_step(func&& A, const Vector & b, const Vector & x, Vector & x_new) {
    // Reuse LU decomposition with SMW
    auto T = A(x);
    Eigen::SparseLU<Eigen::SparseMatrix<double>> Ax_lu;
    Ax_lu.analyzePattern(T); 
    Ax_lu.factorize(T);
    // Solve a bunch of systems
    auto Axinv_b = Ax_lu.solve(b);
    auto Axinv_x = Ax_lu.solve(x);
    // Next step
    x_new = Axinv_b + Ax_lu.solve(x*x.transpose()*(x-Axinv_b)) / (x.norm() + x.dot(Axinv_x) );
}

int main(void) {
    double atol = 1e-13;
    double rtol = 1e-11;
    int max_itr = 100;
    
    // Define a test vector and test rhs and x0 = b
    int n = 8;
    Eigen::SparseMatrix<double> T(n,n);
    T.reserve(3);
    for(int i = 0; i < n; ++i) {
        if(i > 0) T.insert(i,i-1) = 1;
        T.insert(i,i) = 0;
        if(i < n-1) T.insert(i,i+1) = 1;
    }
    
    Eigen::VectorXd b = Eigen::VectorXd::Random(n);
    
    // Define a lambda function implementing A(x)
    // auto = std::function<Eigen::SparseMatrix<double>(const Eigen::VectorXd &)>
    auto A = [&T, n] (const Eigen::VectorXd & x) -> Eigen::SparseMatrix<double> & { double nrm = x.norm();
        for(int i = 0; i < n; ++i) { T.coeffRef(i,i) = 3 + nrm; } return T; };
    
    // Perform convergence study with fixed point iteration
    std::cout << std::endl << "*** Fixed point method ***" << std::endl << std::endl;
    // auto = std::function<Eigen::VectorXd(const Eigen::VectorXd &, Eigen::VectorXd &)>
    auto fix_step = [&A, &b] (const Eigen::VectorXd & x, Eigen::VectorXd & x_new) { fixed_point_step(A, b, x, x_new); };
    
    auto x = b;
    auto x_new = x;
    
    for( int itr = 0;; ) { // Forever until break
        
        // Advance to next step, override x with x_{k+1}
        fix_step(x, x_new);
        
        // Compute residual
        double r = (x - x_new).norm();
        
        std::cout << "[Step " << itr << "] Error: " << r << std::endl;
        
        // Termination conditions
        // If atol reached
        if (r < atol) {
            std::cout << "[CONVERGED] in " << itr << " it. due to atol. err = " << r << " < " << atol << "." << std::endl;
            break;
        }
        // If rtol reached
        if (r < rtol*x_new.norm()) {
            std::cout << "[CONVERGED] in " << itr << " it. due to rtol. err = " << r << " < " << rtol*x_new.norm() << "." << std::endl;
            break;
        }
        // If max it reached
        if (++itr >= max_itr) {
            std::cout << "[NOT CONVERGED] due to MAX it. = " << max_itr << " reached, err = " << r << "." << std::endl;
            break;
        }
        x = x_new;
    }
    
    std::cout << std::endl << "x^*_fix = " << std::endl << x_new << std::endl;
    
    // Perform convergence study with Newton iteration
    std::cout << std::endl << "*** Newton method ***" << std::endl << std::endl;
    
    // auto = std::function<Eigen::VectorXd(const Eigen::VectorXd &, Eigen::VectorXd &)>
    auto newt_step = [&A, &b] (const Eigen::VectorXd & x, Eigen::VectorXd & x_new) { newton_step(A, b, x, x_new); };
    
    x = b;
    
    for( int itr = 0;; ) { // Forever until break
        
        
        // Advance to next step, override x with x_{k+1}
        newt_step(x, x_new);
        
        // Compute residual
        double r = (x - x_new).norm();
        
        std::cout << "[Step " << itr << "] Error: " << r << std::endl;
        
        // Termination conditions
        // If atol reached
        if (r < atol) {
            std::cout << "[CONVERGED] in " << itr << " it. due to atol. err = " << r << " < " << atol << "." << std::endl;
            break;
        }
        // If rtol reached
        if (r < rtol*x_new.norm()) {
            std::cout << "[CONVERGED] in " << itr << " it. due to rtol. err = " << r << " < " << rtol*x_new.norm() << "." << std::endl;
            break;
        }
        // If max it reached
        if (++itr >= max_itr) {
            std::cout << "[NOT CONVERGED] due to MAX it. = " << max_itr << " reached, err = " << r << "." << std::endl;
            break;
        }
        x = x_new;
    }
    
    std::cout << std::endl << "x^*_newt = " << std::endl << x_new << std::endl;
}
