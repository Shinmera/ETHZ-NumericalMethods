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

// Compute the best approximation of the function f with Chebyshev polynomials. alpha is the output vector of coefficients.
template <typename Function>
void bestpolchebnodes(const Function &f, Eigen::VectorXd &alpha) {
    int n=alpha.size()-1;
    Eigen::VectorXd fn(n+1);
    for (int k=0; k<n+1; k++) {
        double temp=cos(M_PI*(2*k+1)/2/(n+1));
        fn(k)=f(temp);
    }
    
    vector<double> V;
    Eigen::MatrixXd scal(n+1,n+1);
    for (int j=0; j<n+1; j++) {
        V=chebpolmult(n,cos(M_PI*(2*j+1)/2/(n+1)));
        for (int k=0; k<n+1; k++) scal(j,k)=V[k];
    }
    
    for (int k=0; k<n+1; k++) {
        alpha(k)=0;
        for (int j=0; j<n+1; j++) {
            alpha(k)+=2*fn(j)*scal(j,k)/(n+1);
        }
    }
        alpha(0)=alpha(0)/2;
}

// Test the implementation.
int main(){
    auto f = [] (double & x) {return 1/(pow(5*x,2)+1);};
    int n=20;
    Eigen::VectorXd alpha(n+1);
    bestpolchebnodes(f, alpha);
    
    //Compute the error
    Eigen::VectorXd X = Eigen::VectorXd::LinSpaced(1e6,-1,1);
    auto qn = [&alpha,&n] (double & x) {
        double temp;
        vector<double> V=chebpolmult(n,x);
        for (int k=0; k<n+1; k++) temp+=alpha(k)*V[k];
        return temp;
    };
    double err_max=abs(f(X(0))-qn(X(0)));
    for (int i=1; i<1e6; i++) err_max=std::max(err_max,abs(f(X(i))-qn(X(i))));
    cout<<"Error: "<< err_max <<endl;
}