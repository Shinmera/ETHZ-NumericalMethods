#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#include <cmath>

#include <iomanip>
#include <iostream>

#define PI  M_PI
#define PI_HALF  M_PI_2

using vector = Eigen::VectorXd;

//! \brief Golub-Welsh implementation 5.3.35
//! \param[in] n number of Gauss nodes
//! \param[out] w weights
//! \param[out] x nodes for interval [-1,1]
void golubwelsh(int n, vector & w, vector & x) {
    w.resize(n);
    x.resize(n);
    if( n == 0 ) {
        x(0) = 0;
        w(0) = 2;
    } else {
        vector b(n-1);
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n,n);

        for(int i = 1; i < n; ++i) {
            double d = (i) / sqrt(4. * i * i - 1.);
            J(i,i-1) = d;
            J(i-1,i) = d;
        }

        Eigen::EigenSolver<Eigen::MatrixXd> eig(J);

        x = eig.eigenvalues().real();
        w = 2 * eig.eigenvectors().real().topRows<1>().cwiseProduct(eig.eigenvectors().real().topRows<1>());
    }
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
    double I = 0;
    for(int i = 0; i < w.size(); ++i) {
        I += f( (x(i) + 1) * (b - a) / 2 + a ) * w(i);
    }
    return I * (b - a) / 2.;
}

//! \brief Compute \int_{-infty}^\infty f(x) dx using transformation x = cot(t)
//! \tparam func template type for function handle f (e.g. lambda func.)
//! \param[in] n number of Gauss points
//! \param[in] f integrand
//! \return Approximation of integral \int_{-infty}^\infty f(x) dx
template <class func>
double quadinf(int n, func&& f) {
    vector w, x;
    golubwelsh(n, w, x);
    // NOTE: no function cot available in c++, need to resort to trigonometric identities
    //! Both below are valid, the first computes two trigonometric functions
    auto ftilde = [&f] (double x) { return f(cos(x)/sin(x)) / pow(sin(x),2); };
//     auto ftilde = [&f] (double x) { double cot = tan(PI_HALF - x); return f(cot) * (1. + pow(cot,2)); };
    return quad(ftilde, w, x, 0, PI);
}

int main() {
    // Number of max Gauss pts.
    double N = 100;

    // Integrand and exact integral
    auto f = [] (double t) { return exp(-pow((t-1),2)); };
    double I_ex = sqrt(PI);
    
    // NOTE: We observe exponential convergence
    int sep = 12;
    std::cout << std::setw(sep) << "Nodes" << std::setw(sep) << "Quadrature" << std::setw(sep) << "Exact" << std::setw(sep) << "Error" << std::endl;
    for(int n = 1; n < N; ++n) {
        vector w, x;
        double QS = quadinf(n, f);
//         std::cout << std::setw(sep) << n << std::setw(sep) << QS << std::setw(sep) << I_ex << std::setw(sep) << std::abs(QS - I_ex) << std::endl;
        std::cout << std::setw(sep) << " " << std::abs(QS - I_ex) << std::endl;
    }
}
