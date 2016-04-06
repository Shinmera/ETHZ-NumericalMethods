#include "taylorintegrator.hpp"

//! \file taylorprey.cpp Solution for Problem 3d, solving prey/predator model with taylor expansion method

int main() {

    // Dimension of state space
    unsigned int d = 2;
    
    // Final time for model
    double T = 10.;
    
    // Initial value for model
    Eigen::VectorXd y0(d);
    y0 << 100, 5;
    
    // Array of number of steps (for convergence study)
    std::vector<unsigned int> N = {128, 256, 512, 1024, 2048, 4096, 8192, 16384};
    
    // Exact value y(10) at final time T = 10 (approximated)
    Eigen::VectorXd yex(d);
    yex << 0.319465882659820, 9.730809352326228;
    
    // Coefficients and handle for prey/predator model
    double alpha1 = 3.;
    double alpha2 = 2.;
    double beta1 = 0.1;
    double beta2 = 0.1;
    // TODO: implement functor f for rhs of y(t)' = f(y(t))
    
    // Constructor
    TaylorIntegrator<Eigen::VectorXd> tint;
    
    // Start convergence study
    std::cout  << std::setw(15) << "N"  << std::setw(15) << "error" << std::setw(15) << "rate" << std::endl;
    // TODO: tabulate error for each N
}
