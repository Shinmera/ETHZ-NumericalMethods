#pragma once

#include <vector>
#include <cassert>

#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

//! \file tylorintegrator.hpp Solution for Problem 3, implementing TaylorIntegrator class

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
        // Iniz step size
        double h = T / N;
        
        // Will contain all steps, reserve memory for efficiency
        std::vector<State> res;
        res.reserve(N);
        
        // Store initial data
        res.push_back(y0);
        
        // Initialize some memory to store temporary values
        State ytemp1 = y0;
        State ytemp2 = y0;
        // Pointers to swap previous value
        State * yold = &ytemp1;
        State * ynew = &ytemp2;
        
        // Loop over all fixed steps
        for(unsigned int k = 0; k < N; ++k) {
            // Compute, save and swap next step
            step(odefun, h, *yold, *ynew);
            res.push_back(*ynew);
            std::swap(yold, ynew);
        }
        
        return res;
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
    void step(const Function &odefun, double h, const State & y0, State & y1) const {
        // Compute values for Taylor expansion, including Jacobian and Hessian matrix
        auto fy0 =  odefun.f(y0);
        auto dfy0fy0 =  odefun.df(y0, fy0);
        auto df2y0fy0 =  odefun.df(y0, dfy0fy0);
        auto d2fy0fy0 = odefun.d2f(y0, fy0);
        
        // Plug values into Taylor expansion for next step
        y1 = y0 + fy0*h + dfy0fy0*h*h/2. + (df2y0fy0 + d2fy0fy0)*h*h*h/6.;
    }
    
};
