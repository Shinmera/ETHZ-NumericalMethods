#include "ode45.hpp"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

// Compute the maps Phi and W at time T, for initial data given by u0 and v0.
pair<Vector2d,Matrix2d> PhiAndW(double u0, double v0, double T) {
    
    auto f = [] (const VectorXd & w) {
        Eigen::VectorXd temp(6);
        temp(0) = (2. - w(1))*w(0);
        temp(1) = (w(0) - 1.)*w(1);
        temp(2) = (2. - w(1))*w(2) - w(0)*w(3);
        temp(3) = w(1)*w(2) + (w(0) - 1.)*w(3);
        temp(4) = (2. - w(1))*w(4) - w(0)*w(5);
        temp(5) = w(1)*w(4) + (w(0) - 1.)*w(5);
        return temp;
    };
    
    Eigen::VectorXd w0(6);
    w0 << u0, v0, 1., 0, 0, 1.;
    
    ode45<Eigen::VectorXd> O(f);
    O.options.rtol = 1e-14;
    O.options.atol = 1e-12;
    auto sol = O.solve(w0, T);
    VectorXd wT = sol.back().first;

    pair<Vector2d,Matrix2d> PaW;
    PaW.first << wT(0), wT(1);
    PaW.second << wT(2), wT(4), wT(3), wT(5);
    return PaW;
}

// Apply the Newton method to find initial data giving solutions with period equal to 5.
int main(){
    Vector2d y;
    y << 3, 2;
    double T = 5;
    pair<Vector2d,Matrix2d> PaW = PhiAndW(y(0), y(1), T);
    Vector2d F = PaW.first - y;
    Matrix2d DF;

    while (F.norm() > 1e-5) {
        PaW = PhiAndW(y(0), y(1), T);
        F = PaW.first - y;
        DF = PaW.second - MatrixXd::Identity(2,2);
        y = y - DF.lu().solve(F);
    }
    
    cout << "The obtained initial condition is: " << endl << y << endl;
    PaW = PhiAndW(y(0), y(1), 100);
    
    cout << "y(100) = " << endl << PaW.first << endl;
}
