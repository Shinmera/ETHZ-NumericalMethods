#include <Eigen/Dense>
#include <utility>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

#include "general_nonlinear_solver.h" // implements general_nonlinear_solver

//! \brief Implements a single step of the modified newton
//! \tparam argument type of argument to function f: such as double or vector etc...
//! \tparam function type for the function f, likely a lambda function
//! \tparam jacobian type for the jacobian df, likely a lambda function
//! \param[in] x previous value to use in step, also initial guess if needed
//! \param[out] x_next next step x_{k+1} of modified Newton
//! \param[in] f function handle for f(x) = 0
//! \param[in] df function handle for jacobian df of f
template <typename argument, class function, class jacobian>
void mod_newt_step_scalar(const argument & x, argument & x_next, function&& f, jacobian&& df) {
    argument y = x + f(x) / df(x);
    x_next = y - f(y) / df(x);
}


//! \brief Implements a single step of the modified newton
//! \tparam argument type of argument to function f: such as double or vector etc...
//! \tparam function type for the function f, likely a lambda function
//! \tparam jacobian type for the jacobian df, likely a lambda function
//! \param[in] x previous value to use in step, also initial guess if needed
//! \param[out] x_next next step x_{k+1} of modified Newton
//! \param[in] f function handle for f(x) = 0
//! \param[in] df function handle for jacobian df of f
template <typename argument, class function, class jacobian>
void mod_newt_step_system(const argument & x, argument & x_next, function&& f, jacobian&& df) {
    auto lu = df(x).lu();
    // reusing LU decomposition in the latter
    // the constructor performs LU decomposition in case argument = Vector
    argument y = x + lu.solve(f(x));
    x_next = y - lu.solve(f(y));
}

//! \brief Solve a scalar non-linear eq. with the modified Newton
void mod_newt_ord() {
    // Setting up values, functions and jacobian
    const double a = 0.123;
    auto f = [&a] (double x) { return atan(x) - a; }; // function, must see a
    auto df = [] (double x) { return 1. / (x*x+1.); };
    
    const double x_ex = tan(a); // exact solution
    
    // store solution and error at each time step
    std::vector<double> sol;
    std::vector<double> err;
    
    // Compute error and push back to err, used in general_nonlinear solver as breaking condition errf(x) < eps
    auto errf = [&err, x_ex] (double & x) { double e =  std::abs(x - x_ex); err.push_back(e); return e; };
    
    // Perform convergence study with Modified newton for scalar
    std::cout << std::endl << "*** Modified Newton method (scalar) ***" << std::endl << std::endl;
    std::cout << "Exact: " << x_ex << std::endl;
    // Initial guess and compute initial error
    double x_scalar = 5;
    sol.push_back(x_scalar);
    errf(x_scalar);
    // Lambda performing the next step, used to define a proper function handle to be passed to general_nonlinear_solver
    auto newt_scalar_step = [&sol, &f, &df] (double x, double & x_new) { mod_newt_step_scalar(x, x_new, f, df); sol.push_back(x_new); };
    // Actually perform the solution
    general_nonlinear_solver(newt_scalar_step, x_scalar, errf);
    // Print solution (final)
    std::cout << std::endl << "x^*_scalar = " << std::endl << x_scalar << std::endl;

    // Print table of solutions, errors and EOOC
    auto space = std::setw(15);
    std::cout << std::endl << space << "sol." << space << "err." << space << "order" << std::endl;
    for(unsigned i = 0; i < sol.size(); ++ i) {
        std::cout << space << sol.at(i) << space << err.at(i);
        if(i >= 3) std::cout << space << (log(err.at(i)) - log(err.at(i-1))) / (log(err.at(i-1)) - log(err.at(i-2)));
        std::cout << std::endl;
    }
}

//! \brief Solve a system non-linear eq. with the modified Newton
void mod_newt_sys() {
    //// Setup initial values
    using Vector = Eigen::VectorXd;
    using Matrix = Eigen::MatrixXd;
    
    Matrix A = Matrix::Random(4,4);
    Vector c = Vector::Random(4);
    
    A = (A*A.transpose()).eval(); // make sure A is symmetric
    c = c.cwiseAbs(); // make sure c is > 0 uniformly
    
    // Handler for function and jacobian in standard format. Must see matrix A and vector c
    auto F = [&A,&c] (const Vector & x) { Vector tmp =  A*x + c.cwiseProduct(x.array().exp().matrix()).eval(); return tmp; };
    std::function<Matrix(const Vector &)> dF = [&A, &c] (const Vector & x) { Matrix C = A; Vector temp = c.cwiseProduct(x.array().exp().matrix()); C += temp.asDiagonal(); return C; };
    
    // Container for errors
    std::vector<double> err;
    // Define lambda for breaking condition, which also stores the error of the previous step. Must see err
    auto rerr = [&err] (Vector & x) { if(err.size() > 0) return err.back(); else return 1.; };
    
    std::cout << std::endl << "*** Modified Newton method (system) ***" << std::endl << std::endl;
    Vector x_system = Vector::Zero(c.size()); // Initial guess
    
    // Refactor the step st. is compatible with general_nonlinear_solver,  must see F; dF and error vector
    auto newt_system_step = [&F, &dF, &err] (const Vector & x, Vector & x_new) { mod_newt_step_system(x, x_new, F, dF); double e = (x - x_new).norm()/ x_new.norm(); err.push_back(e); };
    // Actually performs computations
    general_nonlinear_solver(newt_system_step, x_system, rerr);
    
    // Sanity check
    std::cout << std::endl << "x^*_system = " << std::endl << x_system << std::endl;
    std::cout << std::endl << "F(x^*_system) = " << std::endl << F(x_system) << std::endl;
    
}

int main() {
    // Part 1: solve scalar
    mod_newt_ord();
    
    // Part 2: solve system
    mod_newt_sys();
    
}
