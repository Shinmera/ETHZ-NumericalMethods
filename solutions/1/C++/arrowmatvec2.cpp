#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;

template <class Matrix>

// The auxiliary function Atimesx computes the function A*x in a smart way, using the particular structure of the matrix A.

void Atimesx(const Matrix & d, const Matrix & a, const Matrix & x, Matrix & Ax)
{
    int n=d.size();
    Ax=(d.array()*x.array()).matrix();
    VectorXd Axcut=Ax.head(n-1);
    VectorXd acut = a.head(n-1);
    VectorXd xcut = x.head(n-1);
    
    Ax << Axcut + x(n-1)*acut, Ax(n-1)+ acut.transpose()*xcut;
}

// We compute A*A*x by using the function Atimesx twice with 5 dimensional random vectors.

int main(void)
{
    VectorXd a=VectorXd::Random(5);
    VectorXd d=VectorXd::Random(5);
    VectorXd x=VectorXd::Random(5);
    VectorXd Ax(5);
    
    Atimesx(d,a,x,Ax);
    VectorXd AAx(5);
    Atimesx(d,a,Ax,AAx);
    std::cout << "A*A*x = " << AAx << std::endl;
}
