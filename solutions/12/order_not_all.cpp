#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <iomanip>

#include "errors.hpp"

using namespace std;
using namespace Eigen;

int main() {
    // Construct data for Butcher schemes
    MatrixXd A1 = MatrixXd::Zero(1,1);
    VectorXd b1(1);
    b1 << 1;
    
    MatrixXd A2 = MatrixXd::Zero(2,2);
    A2(1,0) = 1;
    VectorXd b2(2);
    b2 << .5, .5;
        
    MatrixXd A3 = MatrixXd::Zero(3,3);
    A3(1,0) = .5;
    A3(2,0) = -1;
    A3(2,1) = 2;
    VectorXd b3(3);
    b3 << 1./6, 2./3, 1./6;
    
    MatrixXd A4 = MatrixXd::Zero(4,4);
    A4(1,0) = .5;
    A4(2,1) = .5;
    A4(3,2) = 1;
    VectorXd b4(4);
    b4 << 1./6, 1./3, 1./3, 1./6;
    
    // First ODE
    cout << endl << "1. ODE y' = (1-y)y, y(0)=.5" << endl << endl;
    double T = 0.1;
    auto f = [] (VectorXd y) {VectorXd fy(1); fy << (1.-y(0))*y(0); return fy;};
    VectorXd y0(1); y0 << .5;
    
    cout << "Explicit Euler" << endl << endl;
    errors(f, T, y0, A1, b1);
    cout << "Trapezoidal rule" << endl << endl;
    errors(f, T, y0, A2, b2);
    cout << "RK order 3" << endl << endl;
    errors(f, T, y0, A3, b3);
    cout << "Classical RK order 4" << endl << endl;
    errors(f, T, y0, A4, b4);
    
    // Second ODE
    cout << endl << "2. ODE y' = |1.1 - y| + 1, y(0)=1" << endl << endl;
    auto f2 = [] (VectorXd y) {VectorXd fy(1); fy << abs(1.1-y(0))+1.; return fy;};
    y0 << 1;
    
    cout << "Explicit Euler" << endl << endl;
    errors(f2, T, y0, A1, b1);
    cout << "Trapezoidal rule" << endl << endl;
    errors(f2, T, y0, A2, b2);
    cout << "RK order 3" << endl << endl;
    errors(f2, T,  y0, A3, b3);
    cout << "Classical RK order 4" << endl << endl;
    errors(f2, T, y0, A4, b4);
}