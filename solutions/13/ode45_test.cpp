#include "ode45.hpp"

#include <iostream>

#include<Eigen/Dense>
#include <Eigen/Sparse>

// Comment to disable test
#define MAKE_TEST1
#define MAKE_TEST2
#define MAKE_TEST3
#define MAKE_TEST4

void test1() {
#ifdef MAKE_TEST1
    std::cout << "Prey/Predator model test:" << std::endl;
    
    // Prey/Predator model
    const double alpha1 = 3;
    const double alpha2 = 2;
    const double beta1 = 0.1;
    const double beta2 = 0.1;
    auto f = [&alpha1, &alpha2, &beta1, &beta2] (const Eigen::VectorXd & y) {
        Eigen::VectorXd temp = y;
        temp(0) *= (alpha1 - beta1*y(1));
        temp(1) *= (beta2*y(0) - alpha2);
        return temp;
    };
    
//     Eigen::VectorXd y0(2);
    // or
    Eigen::Vector2d y0;
    y0 << 100, 5;
    
    const double T = 10;
    
    // Basic usage:
    ode45<Eigen::Vector2d> O(f);
//     O.options.do_statistics = true;
    auto sol = O.solve(y0, T);
    
    // Print info
//     O.print();
    
    // Print some info
    std::cout << "T = " << sol.back().second << std::endl;
    std::cout << "y(T) = " << std::endl << sol.back().first << std::endl;
#endif
}

void test2() {
#ifdef MAKE_TEST2
    std::cout << "Fundamental types test and validation:" << std::endl;
    
    // Test class with fundamental types
    auto f = [] (double y) { return 1 / y; };
    
    double y0 = 0.2;
    // Large step size
    const double T = 10000;
    
    // Quick syntax
    ode45<double> O(f);
//     O.options.do_statistics = true;
    auto sol = O.solve(y0, T);
    
    // Print info
//     O.print();
    
    // Print some info
    std::cout << "T = " << sol.back().second << std::endl;
    std::cout << "y(T) = " << std::endl << sol.back().first << std::endl;
    auto y_ex = [] (double t) { return std::sqrt(2*t+0.04); };
    std::cout << "y_ex(T) = " << std::endl << y_ex(T) << std::endl;
#endif
}

void test3() {
#ifdef MAKE_TEST3
    std::cout << "Multidimensional test:" << std::endl;
    
    // Construct data for the IVP
    double T = 1;
    // Many dimensions
    int n = 5;
    
    // Multidimensional rhs
    Eigen::VectorXd y0(2*n);
    for(int i = 0; i < n; ++i) {
        y0(i)=(i+1.)/n;
        y0(i+n)=-1;
    }
    
    // Multidimensional rhs
    auto f = [n] (Eigen::VectorXd y) {
        Eigen::VectorXd fy(2*n);
        
        Eigen::VectorXd g(n);
        g(0) = y(0)*(y(1)+y(0));
        g(n-1) = y(n-1)*(y(n-1)+y(n-2));
        for(int i = 1; i < n-1; ++i) {
            g(i) = y(i)*(y(i-1)+y(i+1));
        }
        
        Eigen::SparseMatrix<double> C(n,n);
        C.reserve(3);
        for(int i = 0; i < n; ++i) {
            C.insert(i,i) = 2;
            if(i < n-1) C.insert(i,i+1) = -1;
            if(i >= 1)  C.insert(i,i-1) = -1;
        }
        C.makeCompressed();
        fy.head(n) = y.head(n);
        
        Eigen::SparseLU< Eigen::SparseMatrix<double> >  solver;
        solver.analyzePattern(C);
        solver.compute(C);
        fy.tail(n) = solver.solve(g);
        return fy;
    };
    
    // Constructor:
    ode45<Eigen::VectorXd> O(f);
    
    // Setup options
    O.options.do_statistics = true;
    
    // Solve
    auto sol = O.solve(y0, T);
    
    // Print info
    O.print();

    std::cout << "T = " << sol.back().second << std::endl;
    std::cout << "y(T) = " << std::endl << sol.back().first << std::endl;
#endif
}

void test4() {
#ifdef MAKE_TEST4
    std::cout << "Stiff ode test:" << std::endl;
    
    // Test class with fundamental types
    auto f = [] (double y) { return y*(y - y*y); };
    
    // cf. http://ch.mathworks.com/company/newsletters/articles/stiff-differential-equations.html
    double y0 = 0.0001; // try 0.01
    // Large step size
    const double T = 2 / y0;
    
    // Quick syntax
    ode45<double> O(f);
    O.options.do_statistics = true;
    auto sol = O.solve(y0, T);
    
    // Print info
    O.print();
    
    // Print some info
    std::cout << "T = " << sol.back().second << std::endl;
    std::cout << "y(T) = " << std::endl << sol.back().first << std::endl;
    
    
//     for(auto v: sol) {
//         std::cout << v.first << std::endl;
//     }
#endif
}

int main() {
    
    // Basic prey/predator test
    test1();
    
    // Test double interface
    test2();
    
    // A bit more involved example, also testing options and statistics
    test3();
    
    // A bit more involved example, also testing options and statistics
    test4();
    
    return 0;
}
