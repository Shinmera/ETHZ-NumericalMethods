#pragma once

#include <vector>
#include <cassert>

#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

//! \file rkintegrator.hpp Solution for Problem 1, implementing RkIntegrator class

//! \brief Implements a Runge-Kutta explicit solver for a given Butcher tableau for autonomous ODEs
//! \tparam State a type representing the space in which the solution lies, e.g. R^d, represented by e.g. Eigen::VectorXd.
template <class State>
class RKIntegrator {
public:
    //! \brief Constructor for the RK method.
    //! Performs size checks and copies A and b into internal storage
    //! \param[in] A matrix containing coefficents of Butcher tableau, must be (strictly) lower triangular (no check)
    //! \param[in] b vector containing coefficients of lower part of Butcher tableau
    RKIntegrator(const Eigen::MatrixXd & A, const Eigen::VectorXd & b)
        : A(A), b(b), s(b.size()) {
        assert( A.cols() == A.rows() && "Matrix must be square.");
        assert( A.cols() == b.size() && "Incompatible matrix/vector size.");
    }
    
    //! \brief Perform the solution of the ODE
    //! Solve an autonomous ODE y' = f(y), y(0) = y0, using a RK scheme given in the Butcher tableau provided in the
    //! constructor. Performs N equidistant steps upto time T with initial data y0
    //! \tparam Function type for function implementing the rhs function. Must have State operator()(State x)
    //! \param[in] f function handle for rhs in y' = f(y), e.g. implemented using lambda funciton
    //! \param[in] T final time T
    //! \param[in] y0 initial data y(0) = y0 for y' = f(y)
    //! \param[in] N number of steps to perform. Step size is h = T / N. Steps are equidistant.
    //! \return vector containing all steps y^n (for each n) including initial and final value
    template <class Function>
    std::vector<State> solve(const Function &f, double T, const State & y0, unsigned int N) const {
        // Iniz step size
        double h = T / N;
        
        // Will contain all steps, reserve memory for efficiency
        std::vector<State> res;
        res.reserve(N+1);
        
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
            step(f, h, *yold, *ynew);
            res.push_back(*ynew);
            std::swap(yold, ynew);
        }
        
        return res;
    }
    
private:
    
    //! \brief Perform a single step of the RK method for the solution of the autonomous ODE
    //! Compute a single explicit RK step y^{n+1} = y_n + \sum ... starting from value y0 and storing next value in y1
    //! \tparam Function type for function implementing the rhs. Must have State operator()(State x)
    //! \param[in] f function handle for ths f, s.t. y' = f(y)
    //! \param[in] h step size
    //! \param[in] y0 initial state 
    //! \param[out] y1 next step y^{n+1} = y^n + ...
    template <class Function>
    void step(const Function &f, double h, const State & y0, State & y1) const {
        // create vector holding next value
        y1 = y0;
        
        // Reserve space for increments
        std::vector<State> k;
        k.reserve(s);
        
        // Loop over the size of RK
        for(unsigned int i = 0; i < s; ++i) {
            // Compute increments and save them to k
            State incr = y0;
            for(unsigned int j = 0; j < i; ++j) {
                incr += h*A(i,j)*k.at(j);
            }
            k.push_back( f( incr ) );
            y1 += h * b(i) * k.back();
        }
    }
    
    //! Matrix A in Butcher scheme
    const Eigen::MatrixXd A;
    //! Vector b in Butcher scheme
    const Eigen::VectorXd b;
    //! Size of Butcher matrix and vector A and b
    unsigned int s;
};
