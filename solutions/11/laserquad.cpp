// Include function "gauleg" for the computation of weights/nodes of Gauss quadrature
#include "gauleg.hpp"

#include <iomanip>
#include <iostream>

#include <cmath>

// We have two possible implementations
// A nice one using lambda functions, which is very short but, maybe, more difficult to understand (set this flag)
// A less nice one, using copy-and-paste of the 1D implementation (comment this flag).
#define    USE_LAMBDA

//! \brief Compute \int_a^b f(x) dx \approx \sum w_i f(x_i) (with scaling of w and x)
//! \tparam func template type for function handle f (e.g. lambda func.)
//! \param[in] a left boundary in [a,b]
//! \param[in] b right boundary in [a,b]
//! \param[in] f integrand
//! \param[in] N number of quadrature points
//! \return Approximation of integral \int_a^b f(x) dx
template <class func>
double evalgaussquad(double a, double b, func&& f, const QuadRule & Q) {
    // Compute the qudrature and scale the interal according to inteval [a,b]
    double I = 0;
    for(int i = 0; i < Q.weights.size(); ++i) {
        I += f( (Q.nodes(i) + 1) * (b - a) / 2 + a ) * Q.weights(i);
    }
    return I * (b - a) / 2.;
}

//! \brief Compute double integral \int_\bigtriangleup f(x,b) dx dy using nested Gauss quadrature
//! \tparam func template type for function handle f (e.g. lambda func.), having operator (double x, double y) -> double
//! \param[in] f integrand, f(x,y) must be defined
//! \param[in] N number of quadrature points (in each direction)
//! \return Approximation of integral \int_\bigtriangleup f(x,b) dx dy
template <class func>
double gaussquadtriangle(func&& f, unsigned int N) {
    // Get nodes/weights for integral over dx and dy
    QuadRule Q = gauleg(N);
    
    ////////////////////////////////////////////////
    //// Two liner, lambda-based implementation USE_LAMBDA == true
    ////////////////////////////////////////////////
#ifdef USE_LAMBDA
    // We define an auxiliary function of x defined as f_y := [&y] (double x) { return f(x,y) }, where we fix y
    // We define the function g, as the function of y defined as g(y) := \int_0^{1-y} f_y(x) dx = \int_0^{1-y} I(x,y) dx
    auto g = [&f, &Q] (double y) { return evalgaussquad(0, 1-y, [&f, &y] (double x) { return f(x,y); }, Q); };
    // We integrate the function g over y from 0 to 1
    return evalgaussquad(0, 1, g, Q);
    
    ////////////////////////////////////////////////
    //// Loop based, copy-and-paste implementation USE_LAMBDA == false
    ////////////////////////////////////////////////
#else
    // Integration over y from 0 to 1 of g(y) := \int_0^{1-y} I(x,y) dx
    double I = 0;
    double a = 0., b  = 1.;
    for(int i = 0; i < Q.weights.size(); ++i) {
        // Find out the y at which we are
        double y = (qr.x(i) + 1) * (b - a) / 2 + a;
        // Define f_y(x) (y is fixed and f_y is a function of x)
        auto f_y = [&f, &y] (double x) { return f(x,y); };
        // Compute g(y) as \int_0^{1-y} I(x,y) dx
        I += evalgaussquad(0, 1-y, f_y, Q) * Q.weights(i);
    }
    // Rescale interval
    return I * (b - a) / 2.;
#endif
}

int main(void) {
    // Parameters
    double alpha = 1.;
    double p = 0, q = 0;
    // Laser beam intensity
    auto I = [alpha, p, q] (double x, double y) { return std::exp(- alpha * ( (x-p)*(x-p) + (y-q)*(y-q) ) ); };
    
    // Max num of Gauss points to use (in each direction)
    unsigned int max_N = 13;
    // "Exact" integral
    double I_ex = 0.366046550000405;
    
    // Observed: exponential convergence (as exepcted)
    std::cout << "N" << std::setw(15) << "I_approx" << std::setw(15) << "error" << std::endl;
    for(unsigned int N = 1; N < max_N; ++N) {
        double I_approx = gaussquadtriangle(I, N);
        std::cout << std::setw(3) << N << std::setw(15) << I_approx << std::setw(15) << std::abs(I_ex - I_approx) << std::endl;
    }
}
