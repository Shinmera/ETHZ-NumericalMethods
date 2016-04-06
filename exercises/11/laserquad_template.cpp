// Include function "gauleg" for the computation of weights/nodes of Gauss quadrature
#include "gauleg.hpp"

#include <iomanip>
#include <iostream>

#include <cmath>

//! \brief Compute \int_a^b f(x) dx \approx \sum w_i f(x_i) (with scaling of w and x)
//! \tparam func template type for function handle f (e.g. lambda func.)
//! \param[in] a left boundary in [a,b]
//! \param[in] b right boundary in [a,b]
//! \param[in] f integrand
//! \param[in] Q Structure containing weights and nodes for a quadrature
//! \return Approximation of integral \int_a^b f(x) dx
template <class func>
double evalgaussquad(double a, double b, func&& f, const QuadRule & Q) {
    // TODO: implement gauss quadrature
}

//! \brief Compute double integral \int_\bigtriangleup f(x,b) dx dy using nested Gauss quadrature
//! \tparam func template type for function handle f (e.g. lambda func.), having operator (double x, double y) -> double
//! \param[in] f integrand, f(x,y) must be defined
//! \param[in] N number of quadrature points (in each direction)
//! \return Approximation of integral \int_\bigtriangleup f(x,b) dx dy
template <class func>
double gaussquadtriangle(func&& f, unsigned int N) {
    // TODO: compute double integral over triangle
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
    
    for(unsigned int N = 1; N < max_N; ++N) {
        // TODO: compute intrgral and error
    }
}
