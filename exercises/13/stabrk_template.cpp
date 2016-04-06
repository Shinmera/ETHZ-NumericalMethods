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
    
    
    // TODO: solve IVP of Problem 2 and plot error vs. num. of steps (use uniform timestepping)
}
