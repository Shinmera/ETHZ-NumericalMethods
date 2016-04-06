#include <Eigen/Dense>
#include <iostream>
#include "timer.h"

using namespace Eigen;
using namespace std;

// WARNING: Global! Timer variable
timer<> tm_slow, tm_fast;

//! \brief Compute lmin from the vector d, naive implementation
//! \param[in] d VectorXd of size $n$
//! \param[in] double tol, tolerance
//! \param[out] double lmin

void rankoneinvit(const VectorXd & d, const double & tol, double & lmin)
{
    VectorXd ev=d;
    lmin=0;
    double lnew=d.cwiseAbs().minCoeff();
    while (abs(lnew-lmin)>tol*lmin) {
        tm_slow.start();
        lmin=lnew;
        MatrixXd M=d.asDiagonal();
        M+=ev*ev.transpose();
        ev = M.lu().solve(ev);
        ev.normalize();
        lnew=ev.transpose()*M*ev;
        tm_slow.stop();
    }
    lmin=lnew;
}

//! \brief Compute lmin from the vector d, optimized implementation
//! \param[in] d VectorXd of size $n$
//! \param[in] double tol, tolerance
//! \param[out] double lmin

    
void rankoneinvit_fast(const VectorXd & d, const double & tol, double & lmin)
{
    VectorXd ev=d;
    lmin=0;
    double lnew=d.cwiseAbs().minCoeff();
    
    VectorXd dinv=(1/d.array()).matrix();
    while (abs(lnew-lmin)>tol*lmin) {
        tm_fast.start();
        lmin=lnew;
        VectorXd ev0=ev;
        
        VectorXd Aib=dinv.cwiseProduct(ev); // In these three lines we solve the linear system with the
        double temp=ev.transpose()*Aib;     // Sherman-Morrison-Woodbury formula, in the case of rank-1
        ev=Aib*(1-temp/(1+temp));    // perturbations
        
        
        ev.normalize();
        lnew=ev.transpose()*d.cwiseProduct(ev)+pow(ev.transpose()*ev0,2); //better implementation
        tm_fast.stop();
    }
    lmin=lnew;
}



int main(){
    srand((unsigned int) time(0));
    double tol=1e-3;
    double lmin;
    int n=10;
    
    // Check correctedness of the fast version
    VectorXd d=VectorXd::Random(n);
    rankoneinvit(d,tol,lmin);
    cout<<"lmin = "<<lmin<<endl;
    rankoneinvit_fast(d,tol,lmin);
    cout<<"lmin = "<<lmin<<endl;
    
    //Compare runtimes
    unsigned int repeats = 3;
    
    for(unsigned int p = 2; p <= 9; p++) {
        tm_slow.reset();
        tm_fast.reset();
        unsigned int n = pow(2,p);
        
        for(unsigned int r = 0; r < repeats; ++r) {
            
//             d = VectorXd::Random(n);
            d = VectorXd::LinSpaced(n,1,2);
            rankoneinvit(d,tol,lmin);
            
            rankoneinvit_fast(d,tol,lmin);
        }
        
        cout << "The slow method took:    " <<  tm_slow.min().count() / 1000000. << " ms for n = " <<n<<endl;
        cout << "The fast method took:    " <<  tm_fast.min().count() / 1000000. << " ms for n = " <<n<<endl;
    }
}
