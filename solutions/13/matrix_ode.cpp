#include "ode45.hpp"

#include <iostream>
#include <iomanip>

#include <Eigen/Dense>
#include <Eigen/QR>

//! \file stabrk.cpp Solution for Problem 1, PS13, involving ode45 and matrix ODEs

//! \brief Solve matrix IVP Y' = -(Y-Y')*Y using ode45 up to time T
//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
//! \return Matrix of solution of IVP at t = T
Eigen::MatrixXd matode(const Eigen::MatrixXd & Y0, double T) {

    auto F = [] (const Eigen::MatrixXd & M) { return -(M  - M.transpose())*M; };
    ode45<Eigen::MatrixXd> O(F);
    
    // Set tolerances
    O.options.atol = 10e-10;
    O.options.rtol = 10e-8;
    
    // Return only matrix at T, (solution is vector of pairs (y(t_k), t_k) for each step k
    return O.solve(Y0, T).back().first;
}

//! \brief Find if invariant is preserved after evolution with matode
//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
//! \return true if invariant was preserved (up to round-off), i.e. if norm was less than 10*eps
bool checkinvariant(const Eigen::MatrixXd & M, double T) {
    Eigen::MatrixXd N(3,3);
    
    N = matode(M, T);
    
    if( (N.transpose()*N-M.transpose()*M).norm() < 10 * std::numeric_limits<double>::epsilon() ) {
        return true;
    } else {
        return false;
    }
}

//! \brief Implement ONE step of explicit Euler applied to Y0, of ODE Y' = A*Y
//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
//! \return next step
Eigen::MatrixXd expeulstep(const Eigen::MatrixXd & A, const Eigen::MatrixXd & Y0, double h) {
    return Y0 + h*A*Y0;
}

//! \brief Implement ONE step of implicit Euler applied to Y0, of ODE Y' = A*Y
//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
//! \return next step
Eigen::MatrixXd impeulstep(const Eigen::MatrixXd & A, const Eigen::MatrixXd & Y0, double h) {
    return (Eigen::MatrixXd::Identity(3,3) - h*A).partialPivLu().solve(Y0);
}

//! \brief Implement ONE step of implicit midpoint ruler applied to Y0, of ODE Y' = A*Y
//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
//! \return next step
Eigen::MatrixXd impstep(const Eigen::MatrixXd & A, const Eigen::MatrixXd & Y0, double h) {
    return (Eigen::MatrixXd::Identity(3,3) - h*0.5*A).partialPivLu().solve(Y0+h*0.5*A*Y0);
}

int main() {
    
    double T = 1;
    unsigned int n = 3;
    
    Eigen::MatrixXd M(n,n);
    M << 8,1,6,3,5,7,4,9,2;
    
    std::cout << "SUBTASK 1. c)" << std::endl;
    // Test preservation of orthogonality
    
    // Build Q
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(M.rows(), M.cols());
    qr.compute(M);
    Eigen::MatrixXd Q = qr.householderQ();
    
    // Build A
    Eigen::MatrixXd A(n,n);
    A << 0, 1, 1, -1, 0, 1, -1, -1, 0;
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n,n);
    
    // Norm of Y'Y-I for 20 steps
    Eigen::MatrixXd Mexpeul = Q, Mimpeul = Q, Mimp = Q;
    double h = 0.01;
    std::vector<int> sep = {8,15,15,15};
    std::cout << "Evolution of norm(Y_k'*Y_k - I) for three methods:" << std::endl;
    std::cout   << std::setw(sep[0]) << "step"
                << std::setw(sep[1]) << "exp. Eul"
                << std::setw(sep[2]) << "imp. Eul"
                << std::setw(sep[3]) << "IMP"
                << std::endl;
    std::cout   << std::setw(sep[0]) << "-1"
                << std::setw(sep[1]) << (Mexpeul.transpose()*Mexpeul - I).norm()
                << std::setw(sep[2]) << (Mimpeul.transpose()*Mimpeul - I).norm()
                << std::setw(sep[3]) << (Mimp.transpose()*Mimp - I).norm()
                << std::endl;
    for(unsigned int j = 0; j < 20; ++j) {
        Mexpeul = expeulstep(A, Mexpeul, h);
        Mimpeul = impeulstep(A, Mimpeul, h);
        Mimp = impstep(A, Mimp, h);
        
        std::cout   << std::setw(sep[0]) << j
                    << std::setw(sep[1]) << (Mexpeul.transpose()*Mexpeul - I).norm()
                    << std::setw(sep[2]) << (Mimpeul.transpose()*Mimpeul - I).norm()
                    << std::setw(sep[3]) << (Mimp.transpose()*Mimp - I).norm()
                    << std::endl;
    }
    
    std::cout << "SUBTASK 1. d)" << std::endl;
    // Test implementation of ode45
    
    std::cout << "M = " << std::endl << M << std::endl;
    Eigen::MatrixXd  N = matode(M, T);
    std::cout << "N = " << std::endl << N << std::endl;
    
    std::cout << "SUBTASK 1. g)" << std::endl;
    // Test whether invariant was preserved or not
    
    bool is_invariant = checkinvariant(N, T);
    
    if( is_invariant ) {
        std::cout << "Invariant was preserved." << std::endl;
    } else {
        std::cout << "Invariant was NOT preserved." << std::endl;
    }
    
    
    return 0;
}
