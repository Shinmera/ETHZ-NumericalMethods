#include <Eigen/Dense>

#include <iostream>
#include <iomanip>

#include <cmath>

//! \brief One step of autonomous IVP y'' + y' + y = 0, [y(0), y'(0)] = z0 using SDIRK method
//! Use SDIRK method for first order ode z' = f(z). Steps of size h.
//! \tparam StateType type of solution space y and initial data y0
//! \param[in] z0 initial data z(0)
//! \param[in] h size of the step
//! \param[in] gamma parameter
//! \return next step z1
template <class StateType>
StateType sdirkStep(const StateType & z0, double h, double gamma) {
    // TODO: implement one step of SDIRK method
}

//! \brief Solve autonomous IVP y'' + y' + y = 0, [y(0), y'(0)] = z0 using SDIRK method
//! Use SDIRK method for first order ode z' = f(z), with N equidistant steps
//! \tparam StateType type of solution space z = [y,y']! and initial data z0 = [y(0), y'(0)]
//! \param[in] z0 initial data z(0)
//! \param[in] N number of equidistant steps
//! \param[in] T final time of simulation
//! \param[in] gamma parameter
//! \return vector containing each step of z_k (y and y')
template <class StateType>
std::vector<StateType> sdirkSolve(const StateType & z0, unsigned int N, double T, double gamma) {
    // TODO: implement solution with SDIRK method
}

int main() {
    // Initial data z0 = [y(0), y'(0)]
    Eigen::Vector2d z0;
    z0 << 1,0;
    // Final time
    const double T = 10;
    // Parameter 
    const double gamma = (3.+std::sqrt(3.)) / 6.;
    // Mesh sizes
    std::vector<int> N = {20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240};
    // Exact solution (only y(t)) given z0 = [y(0), y'(0)] and t
    auto yex = [&z0] (double t) {
        return 1./3.*std::exp(-t/2.) * ( 3.*z0(0) * std::cos( std::sqrt(3.)*t/2. ) + 
        std::sqrt(3.)*z0(0) * std::sin( std::sqrt(3.)*t/2. ) + 
        2.*std::sqrt(3.)*z0(1) * std::sin( std::sqrt(3.)*t/2. ) );
    };
    
    // TODO: compute solution and output error/approximate order
}
