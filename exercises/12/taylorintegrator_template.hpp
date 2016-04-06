#pragma once

#include <vector>
#include <cassert>

#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

//! \file tylorintegrator.hpp Solution for Problem 3c, implementing TaylorIntegrator class

//! \brief Implements an autonomous ODE integrator based on Taylor expansion
//! \tparam State a type representing the space in which the solution lies, e.g. R^d, represented by e.g. Eigen::VectorXd.
template <class State>
class TaylorIntegrator {
public:
    //! \brief Perform the solution of the ODE
    //! Solve an autonomous ODE y' = f(y), y(0) = y0, using a Taylor expansion method
    //! constructor. Performs N equidistant steps upto time T with initial data y0
    //! \tparam Function type for function implementing the rhs function (and its derivatives).
    //! \param[in] odefun function handle for rhs f and its derivatives
    //! \param[in] T final time T
    //! \param[in] y0 initial data y(0) = y0 for y' = f(y)
    //! \param[in] N number of steps to perform. Step size is h = T / N. Steps are equidistant.
    //! \return vector containing all steps y^n (for each n) including initial and final value
    template <class Function>
    std::vector<State> solve(const Function &odefun, double T, const State & y0, unsigned int N) const {
        // TODO: solve the autonomous ODE using a Taylor expansion method and suitable call to step
    }
    
private:
    
    //! \brief Perform a single step of the Taylor expansion for the solution of the autonomous ODE
    //! Compute a single explicit step y^{n+1} = y_n + \sum ... starting from value y0 and storing next value in y1
    //! \tparam Function type for function implementing the rhs and its derivatives.
    //! \param[in] odefun function handle for rhs f and the derivatives
    //! \param[in] h step size
    //! \param[in] y0 initial state 
    //! \param[out] y1 next step y^{n+1} = y^n + ...
    template <class Function>
    void step(const Function &odefun, double h,
              const State & y0, State & y1 /* TODO: optional: modify step signature */ ) const {
        // TODO: implement a single step of the Taylor expansion method using provided odefunction
    }
    
};
