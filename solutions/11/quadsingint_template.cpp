 // Include function "gauleg" for the computation of weights/nodes of Gauss quadrature
 #include "gauleg.hpp"
 
 #include <iomanip>
 #include <iostream>
 
 #include <cmath>
 
 //! \brief Compute intrgral \int_{-1}^1 \sqrt(1-t^2) f(t) dt using tranfromation
 //! \param[in] n number of Gauss nodes (evaluate f at 2*n points)
 //! \return value of integral
 template <class func>
 double quadsingint(func&& f, unsigned int n) {
     // TODO: implement quadrature for f(t) * sqrt(1 - t^2) using special transformation and 2*n nodes
 }
 
 int main(void) {
     // Lambda for function f
     auto f = [] (double t) { return 1. / (1. + std::exp(2*t) ); };
     
     // Max num of Gauss points to use (in each direction)
     unsigned int max_N = 25;
     // "Exact" integral
     double I_ex = M_PI_4;
     
     // TODO: estimate q of exponential convergence
     std::cout << "N" << std::setw(15) << "I_approx" << std::setw(15) << "error" << std::setw(15);
     double errold, errnew;
     for(unsigned int N = 1; N < max_N; ++N) {
         double I_approx = quadsingint(f, N);
         
         std::cout << std::setw(3) << N << std::setw(15) << I_approx << std::setw(15) << errnew << std::endl;
     }
     
 }
