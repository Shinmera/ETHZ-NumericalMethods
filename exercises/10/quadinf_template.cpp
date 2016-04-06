#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#include <cmath>

#include <iomanip>
#include <iostream>

#define PI  M_PI

using vector = Eigen::VectorXd;

//! \brief Golub-Welsh implementation 5.3.35
//! \param[in] n number of Gauss nodes
//! \param[out] w weights
//! \param[out] x nodes for interval [-1,1]
void golubwelsh(int n, vector & w, vector & x) {
    // TODO: implement Golub-Welsh
}

//! \brief Compute \int_a^b f(x) dx \approx \sum w_i f(x_i) (with scaling of w and x)
//! \tparam func template type for function handle f (e.g. lambda func.)
//! \param[in] f integrand
//! \param[in] w weights
//! \param[in] x nodes for interval [-1,1]
//! \param[in] a left boundary in [a,b]
//! \param[in] b right boundary in [a,b]
//! \return Approximation of integral \int_a^b f(x) dx
template <class func>
double quad(func&& f, const vector & w, const vector & x, double a, double b) {
    // TODO: implement generic quadrature
    // WARNING: careful scaling
}

//! \brief Compute \int_{-infty}^\infty f(x) dx using transformation x = cot(t)
//! \tparam func template type for function handle f (e.g. lambda func.)
//! \param[in] n number of Gauss points
//! \param[in] f integrand
//! \return Approximation of integral \int_{-infty}^\infty f(x) dx
template <class func>
double quadinf(int n, func&& f) {
    // TODO: implement tranformation of f and call to quad
}

int main() {
    // TODO: test integration of h with 1 to 100 Gaussian points
}
