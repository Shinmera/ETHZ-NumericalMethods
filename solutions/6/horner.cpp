#include <vector>
#include <complex>
#include <iostream>
#include "timer.h"

using namespace std;

// WARNING: Global! Timer variable
timer<> tm_slow, tm_fast;

//! \brief Evaluate a polynomial and its derivative using Horner scheme
//! \param[in] vector c of size $n$, coefficients of the polynomial p
//! \param[in] double x, where the polynomial has to be evaluated
//! \param[out] pair containing p(x),p'(x)

template <typename CoeffVec>
std::pair<double, double> evaldp ( const CoeffVec & c, const double x )
{
    std::pair<double, double> p;
    double px,dpx;
    int s=c.size();

    px = c[0];
    for (int i=1; i<s; i++) {px = x*px+c[i];}
    
    dpx=(s-1)*c[0];
    for (int i=1; i<s-1; i++) {dpx = x*dpx+(s-i-1)*c[i];}
    
    p.first=px;
    p.second=dpx;
    return p;
}

//! \brief Evaluate a polynomial and its derivative using a naive implementation
//! \param[in] vector c of size $n$, coefficients of the polynomial p
//! \param[in] double x, where the polynomial has to be evaluated
//! \param[out] pair containing p(x),p'(x)

template <typename CoeffVec>
std::pair<double, double> evaldp_naive ( const CoeffVec & c, const double x )
{
    std::pair<double, double> p;
    double px,dpx;
    int n=c.size();
    
    px = c[0]*pow(x,n-1);
    for (int i=1; i<n; i++) {px = px+c[i]*pow(x,n-i-1);}
    
    dpx=(n-1)*c[0]*pow(x,n-2);
    for (int i=1; i<n-1; i++) {dpx = dpx+(n-i-1)*c[i]*pow(x,n-i-2);}
    
    p.first=px;
    p.second=dpx;
    return p;
}

int main(){
    vector<double> c {3, 1, 5, 7, 9};
    double x=.123;
    
    // Check implementations
    std::pair<double, double> p,p_naive;
    p=evaldp(c,x);
    cout<<p.first<<endl<<p.second<<endl;
    
    p_naive=evaldp_naive(c,x);
    cout<<p_naive.first<<endl<<p_naive.second<<endl;
    
    //Compare runtimes
    unsigned int repeats = 10;
    
    for(unsigned int k = 2; k <= 20; k++) {
        tm_slow.reset();
        tm_fast.reset();
        vector<double> c;
        
        int n=pow(2,k);
        for (int i=0; i<n; i++) {c.push_back(i+1);}

            for(unsigned int r = 0; r < repeats; ++r) {
            
             tm_slow.start();
             p_naive=evaldp_naive(c,x);
             tm_slow.stop();
                
             tm_fast.start();
             p=evaldp(c,x);
             tm_fast.stop();
        }
        
        cout << "The slow method took:    " <<  tm_slow.avg().count() / 1000000. << " ms for n = " <<n<<endl;
        cout << "The fast method took:    " <<  tm_fast.avg().count() / 1000000. << " ms for n = " <<n<<endl;
    }
}