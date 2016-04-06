// Include function "gauleg" for the computation of weights/nodes of Gauss quadrature
#include "gauleg.hpp"

#include <iomanip>
#include <iostream>

#include <cmath>

// Choose between substitution 1. s = sqrt(1 \pm t), 2. t = cos(s)
#define METHOD 2

//! \brief Compute intrgral \int_{-1}^1 \sqrt(1-t^2) f(t) dt using tranfromation
//! \param[in] n number of Gauss nodes (evaluate f at 2*n points)
//! \return value of integral
template <class func>
double quadsingint(func&& f, unsigned int n) {
    double I = 0.;
#if METHOD == 1 // s = sqrt(1 \pm t)
    QuadRule Q = gauleg(n);
    
    for(unsigned int i = 0; i < n; ++i) {
        double x = (Q.nodes(i) + 1.) / 2.;
        double w = Q.weights(i) * x * x * std::sqrt(2. - x * x);
        I += w * (f(x*x - 1) + f(-x*x + 1));
    }
    return I;
#elif METHOD == 2 // t = cos(s)
    QuadRule Q = gauleg(2*n);
    
    for(unsigned int i = 0; i < 2*n; ++i) {
        double x = sin(Q.nodes(i) * M_PI_2);
        double w = Q.weights(i) * cos(Q.nodes(i) * M_PI_2) * cos(Q.nodes(i) * M_PI_2);
        I += w * f(x) * M_PI_2;
    }
    return I;
#endif
}

int main(int argc, char ** argv) {
    
    auto f = [] (double t) { return 1. / (2. + std::exp(3*t) ); };
    
    // Max num of Gauss points to use (in each direction)
    unsigned int max_N = 25;
    // "Exact" integral
    double I_ex = 0.483296828976607;
    
    // IF you want to test monomials as f, use this data (uncomment this and comment f above:
//     std::vector<double> ex = { M_PI_2, 0., M_PI_2 / 4, 0, M_PI_2 / 8, 0, M_PI_2 * 5 / 64};
//     assert(argc > 1);
//     int n = std::atoi(argv[1]);
//     std::cout << "Monomial of degree:" << n << std::endl;
//     auto f = [&n] (double t) { return std::pow(t, n); };
//     I_ex = ex[n];
    
    // Observed: exponential convergence (as exepcted)
    std::cout << "N" << std::setw(15) << "I_approx" << std::setw(15) << "error" << std::setw(15) << "q est." << std::endl;
    double errold, errnew;
    for(unsigned int N = 1; N < max_N; ++N) {
        double I_approx = quadsingint(f, N);
        
        errold = errnew;
        errnew = std::abs(I_ex - I_approx);
        std::cout << std::setw(3) << N << std::setw(15) << I_approx << std::setw(15) << errnew;
        
        
        if( N > 1 ) {
            std::cout  << std::setw(15) << errnew / errold;
        }
        std::cout << std::endl;
    }

}
