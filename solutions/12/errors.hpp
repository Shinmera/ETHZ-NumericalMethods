#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <iomanip>

#include "rkintegrator.hpp"

using namespace std;
using namespace Eigen;

// This function approximates the order of convergence of the RK scheme defined by A and b when applied to the first order system y'=f(y), y(0)=y0. We are interested in the error of the solutions at the point T.

template <class Function>
void errors(const Function &f, const double &T, const VectorXd &y0, const MatrixXd &A, const VectorXd &b) {
    
    RKIntegrator<VectorXd> rk(A,b);
    vector<double> error(15);
    vector<double> order(14);
    double sum = 0;
    int count = 0;
    bool test = 1;
    vector<VectorXd> y_exact = rk.solve(f,T,y0,pow(2,15));
    
    for(int k = 0; k < 15; k++) {
        int N = pow(2,k+1);
        vector<VectorXd> y1 = rk.solve(f,T,y0,N);
        
        error[k] = (y1[N]-y_exact[pow(2,15)]).norm();
        cout << left << setw(3) << setfill(' ') << "N = ";
        cout << left << setw(7) << setfill(' ') << N;
        cout << left << setw(8) << setfill(' ') << "Error = ";
        cout << left << setw(13) << setfill(' ') << error[k];
        
        if (error[k]<y0.size()*5e-14) test = 0;
        if (k>0 && test) {
            order[k-1]=log(error[k-1]/error[k])/log(2);
            cout << left << setw(10) << setfill(' ') << "Approximated order = " << order[k-1] <<endl;
            sum += order[k-1];
            count = k;
        }
        else cout << endl;
    }
    cout << "Average approximated order = " << sum / count << endl << endl;
}
