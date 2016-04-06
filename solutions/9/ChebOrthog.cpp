#include <iostream>
#include <math.h>
#include <vector>
#include <Eigen/Dense>

using namespace std;

//Evaluate the Chebyshev polynomials up to order n in x.
vector<double> chebpolmult(const int &n,const double &x)
{
    vector<double> V={1,x};
    for (int k=1; k<n; k++)
        V.push_back(2*x*V[k]-V[k-1]);
    return V;
}

//Check orthogonality of Chebyshev polynomials.
int main(){
    int n=10;
    vector<double> V;
    Eigen::MatrixXd scal(n+1,n+1);
    for (int j=0; j<n+1; j++) {
        V=chebpolmult(n,cos(M_PI*(2*j+1)/2/(n+1)));
        for (int k=0; k<n+1; k++) scal(j,k)=V[k];
    }
    cout<<"Scalar products: "<<endl;
    for (int k=0; k<n+1; k++)
        for (int l=k+1; l<n+1; l++)
            cout<<scal.col(k).dot(scal.col(l))<<endl;
}