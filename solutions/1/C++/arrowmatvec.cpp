#include <Eigen/Dense>
#include <iostream>
#include <ctime>

using namespace Eigen;

template <class Matrix>

// The function arrowmatvec computes the product the desired product directly as A*A*x
void arrowmatvec(const Matrix & d, const Matrix & a, const Matrix & x, Matrix & y)
{
    // Here, we choose a MATLAB style implementation using block construction
    // you can also use loops
    // If you are interested you can compare both implementation and see if and how they differ
    int n=d.size();
    VectorXd dcut= d.head(n-1);
    VectorXd acut = a.head(n-1);
    MatrixXd ddiag=dcut.asDiagonal();
    MatrixXd A(n,n);
    MatrixXd D = dcut.asDiagonal();
    // If you do not create the temporary matrix D, you will get an error: D must be casted to MatrixXd
    A << D, a.head(n-1), acut.transpose(), d(n-1);
    
    y=A*A*x;
}

// We test the function arrowmatvec with 5 dimensional random vectors.

int main(void)
{
//     srand((unsigned int) time(0));
    VectorXd a=VectorXd::Random(5);
    VectorXd d=VectorXd::Random(5);
    VectorXd x=VectorXd::Random(5);
    VectorXd y;
    
    arrowmatvec(d,a,x,y);
    std::cout << "A*A*x = " << y << std::endl;
}
