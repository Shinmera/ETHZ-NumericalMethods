#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <vector>
#include <iomanip>

using namespace Eigen;
using namespace std;

// This C++ function phim gives the function phi used for the construction of the exponential Euler single step method for an autonomous ODE.
MatrixXd  phim(MatrixXd Z) {
    int n = Z.cols();
    assert( n == Z.rows() && "Matrix must be square.");
    MatrixXd C(2*n,2*n);
    C << Z, MatrixXd::Identity(n,n), MatrixXd::Zero(n,2*n);
    return C.exp().block(0,n,n,n);
}

// This function calculates a single step of the exponential Euler method, where y0 is the initial state, f and df are object with evaluation operators representing f and df, and h is the stepsize.
template <class Function, class Function2>
VectorXd ExpEulStep(VectorXd y0, Function f, Function2 df, double h) {
   // TODO
}

// Test the exponential Euler method with the logistic ODE and determine the approximated order of algebraic convergence.
int main() {
    // TODO
}