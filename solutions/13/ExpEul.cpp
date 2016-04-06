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
   return  y0 + h * phim( h * df(y0)) * f(y0);
}

// Test the exponential Euler method with the logistic ODE and determine the approximated order of algebraic convergence.
int main() {
    double T = 1;
    VectorXd y0(1); y0 << 0.1;
    auto f = [] (VectorXd y) {return y(0)*(1-y(0));};
    auto df = [] (VectorXd y) {VectorXd dfy(1); dfy << 1-2*y(0); return dfy;};
    double exactyT = y0(0)/(y0(0)+(1-y0(0))*exp(-T));
    
    vector<double> error(15);
    
    for (int j=0; j < 15; j++) {
        int N = pow(2,j+1);
        double h = T / N;
        VectorXd y = y0;
        for (int k=0; k < N; k++) y = ExpEulStep(y,f,df,h);
        
        error[j] = abs(y(0) - exactyT);
        cout << left << setw(3) << setfill(' ') << "N = ";
        cout << left << setw(7) << setfill(' ') << N;
        cout << left << setw(8) << setfill(' ') << "Error = ";
        cout << left << setw(13) << setfill(' ') << error[j];
        if (j > 0)  cout << left << setw(10) << setfill(' ') << "Approximated order = " << log(error[j-1]/error[j])/log(2) <<endl;
        else cout << endl;
    }
}