#include "ode45.hpp"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

// Compute the maps Phi and W at time T, for initial data given by u0 and v0.
pair<Vector2d,Matrix2d> PhiAndW(double u0, double v0, double T) {
    auto f = [] (const VectorXd & w) {
        
        // TODO: the right hand side of the system of ODEs related to Phi and W.
        
    };
    
    ode45<Eigen::VectorXd> O(f);
    O.options.rtol = 1e-14;
    O.options.atol = 1e-12;

    // TODO
    
}

// Apply the Newton method to find initial data giving solutions with period equal to 5.
int main(){
    
    // TODO
    
}