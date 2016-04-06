#include <iostream>
#include <math.h>
#include <vector>
#include <complex>

using namespace std;

// Function implementing weighted Gauss quadrature
template<typename Function>
double quadU(const Function &f, unsigned int n) {
    double q = 0;
    double w,xi;
    for (int j = 0; j < n; j++) {
        w = M_PI/(n+1)*pow(sin((j+1)*M_PI/(n+1)),2);
        xi = cos((j+1.)/(n+1)*M_PI);
        q += w*f(xi);
    }
    return q;
}


// Test the implementation. The parameter q of the exponential decay is approximated by 0.1-0.2. After the 18th iteration only numerical error is present.
int main(){
    auto f = [] (double & t) {return 1/(2+exp(3*t));};
    double exact = 0.483296828976607;
    vector<double> e(25);
    for (unsigned int n = 0; n < 25; n++) {
        e[n] = abs(quadU(f,n+1)-exact);
        if (n>1)
            cout<<"Error with n="<<n<<": "<<e[n]<<". Approximated q: "<<e[n]/e[n-1]<<endl;
    }
}