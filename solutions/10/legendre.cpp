#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Evaluate the Legendre polynomials and its derivatives in the vector x using the 3-term recursion formulae. The outputs are the matrices Lx and DLx.
void legvals(const VectorXd &x, MatrixXd &Lx, MatrixXd &DLx)
{
    int n = Lx.cols()-1;
    int N = x.size();
    for (int j=0; j<N; j++) {
        Lx(j,0) = 1.;
        Lx(j,1) = x(j);
        DLx(j,0) = 0;
        DLx(j,1) = 1.;
        for (int k=2; k<n+1; k++) {
            Lx(j,k) = (2*k-1.)/k*x(j)*Lx(j,k-1)-(k-1.)/k*Lx(j,k-2);
            DLx(j,k) = (2*k-1.)/k*Lx(j,k-1)+(2*k-1.)/k*x(j)*DLx(j,k-1)-(k-1.)/k*DLx(j,k-2);
        }
    }
}

// Evaluate P_n(x), for a scalar x and integer n.
double Pnx(double x, int n) {
    VectorXd Px(n+1);
    Px(0) = 1.; Px(1) = x;
    for (int k=2; k<n+1; k++)
        Px(k) = (2*k-1.)/k*x*Px(k-1)-(k-1.)/k*Px(k-2);
    return Px(n);
}

// Find the Gauss points using the secant method with regula falsi. The standard secant method may be obtained by commenting out lines 50 and 52.
MatrixXd gaussPts(int n, double rtol=1e-10, double atol=1e-12) {
    MatrixXd zeros(n,n);
    double x0, x1, f0, f1, s;
    for (int k=1; k<n+1; k++) {
        for (int j=1; j<k+1; j++) {
            // Initialise initial guesses.
            if (j==1) x0 = -1.;
            else      x0 = zeros(j-2,k-2);
            if (j==k) x1 = 1.;
            else      x1 = zeros(j-1,k-2);
            
            // Secant method
            f0 = Pnx(x0,k);
            for (int i=0; i<1e4; i++) {
                f1 = Pnx(x1,k);
                s = f1*(x1-x0)/(f1-f0);
                if (Pnx(x1 - s,k)*f1<0) {
                    x0 = x1; f0 = f1;
                }
                x1 = x1 - s;
                if ((abs(s)<max(atol,rtol*min(abs(x0),abs(x1)))))  {
                    zeros(j-1,k-1) = x1;
                    break;
                }
            }
        }
    }
    return zeros;
}

// Test the implementation.
int main(){
    int n = 8;
    MatrixXd zeros = gaussPts(n);
    cout<<"Zeros: "<<endl<< zeros <<endl;
    
    for (int k=1; k<n+1; k++) {
        VectorXd xi = zeros.block(0, k-1, k, 1);
        MatrixXd Lx(k,n+1), DLx(k,n+1);
        legvals(xi, Lx, DLx);
        cout<<"Values of the "<<k<<"-th polynomial in the calculated zeros: "<<endl;
        cout<<Lx.col(k).transpose() <<endl;
    }
}
