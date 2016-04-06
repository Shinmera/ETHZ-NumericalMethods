#include "implicit_rkintegrator.hpp"

//! \file cross.cpp Solution for Problem 2, solving cross-product ODE with implicit RK method


//! \brief Perform the solution of the ODE with the linear implicit mid-point method
//! Solve an autonomous ODE y' = f(y), y(0) = y0, using the linear implicit mid-point method. Performs N equidistant steps upto time T with initial data y0
//! \tparam Function type for function implementing the rhs function. Must have Eigen::VectorXd operator()(Eigen::VectorXd x)
//! \tparam Function2 type for function implementing the Jacobian of f. Must have Eigen::MatrixXd operator()(Eigen::VectorXd x)
//! \param[in] f function handle for rhs in y' = f(y), e.g. implemented using lambda funciton
//! \param[in] Jf function handle for Jf, e.g. implemented using lambda funciton
//! \param[in] T final time T
//! \param[in] y0 initial data y(0) = y0 for y' = f(y)
//! \param[in] N number of steps to perform. Step size is h = T / N. Steps are equidistant.
//! \return vector containing all steps y^n (for each n) including initial and final value
template <class Function, class Function2>
std::vector<Eigen::VectorXd> solve_lin_mid(const Function &f, const Function2 &Jf, double T, const Eigen::VectorXd & y0, unsigned int N)  {
    // TODO
}


int main(void) {
    
    // 1. Implicit mid-point method
    
    // TODO: solve the cros-product ODE with the class implicit_RKIntegrator and the implicit mid-point method.
    
    // 2. Linear mid-point method
    
    // TODO: solve the cros-product ODE with the function solve_lin_mid.