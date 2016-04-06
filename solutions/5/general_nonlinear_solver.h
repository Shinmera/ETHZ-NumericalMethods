#pragma once

//! \param[in] res function implementing the step x_{k+1} = step(x_{k}), signature step(x_{k}, x_{k+1})
//! \param[in,out] x initial data (as input) and final iteration (as output)
//! \param[in] errf function implementing the norm of the error (errf(x)) for breaking condition
//! \param[in] eps tolerance to break iterations when res(x) < eps
//! \param[in] max_itr maximal number of iterations
template <class stp_func, class Vector, class err_func>
bool general_nonlinear_solver(stp_func&& step, Vector & x, err_func&& errf, double eps = 1e-8, int max_itr = 100) {
    // Temporary where to store new step
    auto x_new = x;
    
    for( int itr = 0;; ) { // Forever until break
        
        // Compute residual
        double r = errf(x);
        
        std::cout << "[Step " << itr << "] Error: " << r << std::endl;
        
        // Advance to next step, override x with x_{k+1}
        step(x, x_new);
        
        // Termination conditions
        // If tol reached
        if (r < eps) {
            std::cout << "[CONVERGED] in " << itr << " it. due to err. err = " << r << " < " << eps << "." << std::endl;
            return true;
        }
        // If max it reached
        if (++itr >= max_itr) {
            std::cout << "[NOT CONVERGED] due to MAX it. = " << max_itr << " reached, err = " << r << "." << std::endl;
            return false;
        }
        x = x_new;
    }
}
