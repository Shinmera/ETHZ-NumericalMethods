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
    
    // TODO: implement rosenbrock method
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

    // Function and his Jacobian
    // TODO: implement r.h.s. and Jacobian
        
    // TODO: compute reference solutions, solution, error and output approximate order of convergence
    }
}
                                    
