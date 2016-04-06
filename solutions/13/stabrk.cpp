#include "rkintegrator.hpp"

//! \file stabrk.cpp Solution for Problem 2, PS13, solving prey/predator model with RK-SSM method

int main(void) {
    
    // Implementation of butcher scheme
    unsigned int s = 3;
    Eigen::MatrixXd A(s,s);
    Eigen::VectorXd b(s);
    A << 0,      0,      0,
         1.,     0,      0,
         1./4.,  1./4.,  0;
    b << 1./6.,  1./6.,  2./3.;
    
    
    // Coefficients and handle for prey/predator model
    double alpha1 = 1.;
    double alpha2 = 1.;
    double beta1 = 1.;
    double beta2 = 1.;
    auto f = [&alpha1, &alpha2, &beta1, &beta2] (const Eigen::VectorXd & y) { 
        Eigen::VectorXd temp = y;
        temp(0) *= alpha1 - beta1*y(1);
        temp(1) *= -alpha2 + beta2*y(0);
        return temp;
    };
    
    // Dimension of state space
    unsigned int d = 2;
    
    // Final time for model
    double T = 1.;
    
    // Initial value for model
    Eigen::VectorXd y0(d);
    y0 << 100, 1;
    
    // Array of number of steps (for convergence study)
    std::vector<unsigned int> N = {4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192};
    unsigned int k = 16384;
    
    // Initialize RK with Butcher table
    RKIntegrator<Eigen::VectorXd> RK(A,b);
    auto sol  = RK.solve(f, T, y0, k);
    for(auto v: sol) {
        std::cout << v << std::endl;
    }
    
    // Exact value y(10) at final time T (approximated)
    Eigen::VectorXd yex  = RK.solve(f, T, y0, k).back();
    
    // Start convergence study
    std::cout  << std::setw(15) << "N"  << std::setw(15) << "error" << std::setw(15) << "rate" << std::endl;
    double errold;
    for(unsigned int i = 0; i < N.size(); ++i) {
        auto res = RK.solve(f, T, y0, N[i]);
        double err = (res.back() - yex).norm();
        std::cout  << std::setw(15) << N[i] << std::setw(15) << err;
        if(i > 0) {
            std::cout << std::setw(15) << log2(errold / err);
        }
        errold = err;
        std::cout  << std::endl;
    }
}
