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
    // Iniz step size
    double h = T / N;
    int d = y0.size();
    
    // Will contain all steps, reserve memory for efficiency
    std::vector<Eigen::VectorXd> res;
    res.reserve(N+1);
    
    // Store initial data
    res.push_back(y0);
    
    // Initialize some memory to store temporary values
    Eigen::VectorXd ytemp1 = y0;
    Eigen::VectorXd ytemp2 = y0;
    // Pointers to swap previous value
    Eigen::VectorXd * yold = &ytemp1;
    Eigen::VectorXd * ynew = &ytemp2;
    
    // Loop over all fixed steps
    for(unsigned int k = 0; k < N; ++k) {
        // Compute, save and swap next step
        *ynew = *yold+h*(Eigen::MatrixXd::Identity(3,3)-h/2.*Jf(*yold)).lu().solve(f(*yold));
        res.push_back(*ynew);
        std::swap(yold, ynew);
    }
    return res;
}


int main(void) {
    
    // 1. Implicit mid-point method
    std::cout << "1. Implicit midpoint method" << std::endl << std::endl;
    
    unsigned int s = 1;
    Eigen::MatrixXd A(s,s);
    Eigen::VectorXd b(s);
    A << 1./2.;
    b << 1.;
    
    // Coefficients and handle for ODE
    double c = 1.;
    Eigen::Vector3d a;
    a << 1., 0., 0.;
    auto f = [&a, &c] (const Eigen::Vector3d & y) -> Eigen::Vector3d {
        return a.cross(y) + c*y.cross(a.cross(y));
        };
    
   // auto f = [&a, &c] (const Eigen::Vector3d & y) {
   //     Eigen::Vector3d temp;
   //     temp <<   a(1)*y(2) - a(2)*y(1) + c*(y(1)*(a(0)*y(1) - a(1)*y(0)) + y(2)*(a(0)*y(2) - a(2)*y(0))),
   //     a(2)*y(0) - a(0)*y(2) - c*(y(0)*(a(0)*y(1) - a(1)*y(0)) - y(2)*(a(1)*y(2) - a(2)*y(1))),
   //     a(0)*y(1) - a(1)*y(0) - c*(y(0)*(a(0)*y(2) - a(2)*y(0)) + y(1)*(a(1)*y(2) - a(2)*y(1)));
   //     return temp;
   // };
    
    
    auto Jf = [&a, &c] (const Eigen::Vector3d & y) {
        Eigen::Matrix3d temp;
        temp << -c*(a(1)*y(1) + a(2)*y(2)),           c*(2*a(0)*y(1) - a(1)*y(0)) - a(2),  a(1) + c*(2*a(0)*y(2) - a(2)*y(0)),
                a(2) - c*(a(0)*y(1) - 2*a(1)*y(0)),   -c*(a(0)*y(0) + a(2)*y(2)),          c*(2*a(1)*y(2) - a(2)*y(1)) - a(0),
                - a(1) - c*(a(0)*y(2) - 2*a(2)*y(0)), a(0) - c*(a(1)*y(2) - 2*a(2)*y(1)),  -c*(a(0)*y(0) + a(1)*y(1));
        return temp;
    };
    
    // Initial value and final time
    Eigen::Vector3d y0;
    y0 << 1., 1., 1.;
    double T = 10.;
   
    // Initialize implicit RK with Butcher scheme
    implicit_RKIntegrator RK(A,b);
    int N = 128;
    
    auto res = RK.solve(f, Jf, T, y0, N);
    for (int i=0; i<N+1; i++) {
        std::cout<< "norm(y(" << T*i/N << ")) = " << res[i].norm()<<std::endl;
    }
    
    // 2. Linear mid-point method
    std::cout << std::endl << "2. Linear implicit midpoint method" << std::endl << std::endl;
    
    res = solve_lin_mid(f, Jf, T, y0, N);
    for (int i=0; i<N+1; i++) {
        std::cout<< "norm(y(" << T*i/N << ")) = " << res[i].norm()<<std::endl;
    }
}