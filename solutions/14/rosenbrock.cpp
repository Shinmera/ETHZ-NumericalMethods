#include <Eigen/Dense>

#include <iostream>
#include <iomanip>

#include <cmath>

//! \brief Solve the autonomous IVP y' = f(y), y(0) = y0 using Rosenbrock method
//! Use semi-implicit Rosenbrock method using Jacobian evaluation. Equidistant steps of size T/N.
//! \tparam Func function type for r.h.s. f
//! \tparam DFunc function type for Jacobian df
//! \tparam StateType type of solution space y and initial data y0
//! \param[in] f r.h.s. func f
//! \param[in] df Jacobian df of f
//! \param[in] y0 initial data y(0)
//! \param[in] N number of equidistant steps
//! \param[in] T final time
//! \return vector of y_k for each step k from 0 to N
template <class Func, class DFunc, class StateType>
std::vector<StateType> solveRosenbrock(const Func & f, const DFunc & df,
                                       const StateType & y0, unsigned int N, double T) {
    const double h = T/N;
    const double a = 1. / (std::sqrt(2) + 2.);
    
    // Will contain all time steps
    std::vector<StateType> res(N+1);
    // Push initial data
    res.at(0) = y0;
    
    // Some temporary
    StateType k1, k2;
    Eigen::Matrix2d J, W;
    
    // Main loop: *up to N (performs N steps)*
    for(unsigned int i = 1; i <= N; ++i) {
        StateType & yprev = res.at(i-1);
        
        // Jacobian computation
        J = df(yprev);
        W = Eigen::Matrix2d::Identity() - a*h*J;
        
        // Reuse factorization for each step
        auto W_lu = W.partialPivLu();
        
        // Increments
        k1 = W_lu.solve(f(yprev));
        k2 = W_lu.solve(f(yprev+0.5*h*k1) - a*h*J*k1);
        
        // Push new step
        res.at(i) = yprev + h*k2;
    }
    
    return res;
}


int main() {
    // Final time
    const double T = 10;
    // All mesh sizes
    const std::vector<int> N = {16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
    // Reference mesh size
    const int N_ref = 16384;
    // Initial data
    Eigen::Vector2d y0;
    y0 << 1., 1.;
    // Parameter and useful matrix for f
    const double lambda = 1;
    Eigen::Matrix2d R;
    R << 0., -1., 1., 0.;

    // Function and his Jacobian
    auto f =  [&R, &lambda] (const Eigen::Vector2d & y) { return R*y + lambda*(1. - std::pow(y.norm(),2))*y; };
    auto df = [&lambda] (const Eigen::Vector2d & y) {
        double x = 1 - std::pow(y.norm(), 2);
        Eigen::Matrix2d J;
        J << lambda*x - 2*lambda*y(0)*y(0),
        -1 - 2*lambda*y(1)*y(0),
        1 - 2*lambda*y(1)*y(0),
        lambda*x - 2*lambda*y(1)*y(1);
        return J;
    };

    // Reference solution
    auto solref = solveRosenbrock(f, df, y0, N_ref, T);

    std::cout << std::setw(15) << "n" << std::setw(15) << "maxerr" << std::setw(15) << "rate" << std::endl;
    // Used to compute rate at each step:
    double errold = 0;
    // Main loop: loop over all meshes
    for(unsigned int i = 0; i < N.size(); ++i) {
        int n = N.at(i);
        
        // Get solution
        auto sol = solveRosenbrock(f, df, y0, n, T);
        
        // Compute error
        double maxerr = 0;
        for(unsigned int j = 0; j < sol.size(); ++j) {
            maxerr = std::max(maxerr, (sol.at(j) - solref.at((j*N_ref)/n)).norm());
        }
        
        // Do some I/O
        std::cout << std::setw(15) << n << std::setw(15) << maxerr;
        if(i > 0) std::cout << std::setw(15) << std::log2(errold / maxerr);
        std::cout << std::endl;
        
        // Save old error (for rate)
        errold = maxerr;
    }
}
                                    
