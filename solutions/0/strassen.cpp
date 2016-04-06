#include <Eigen/Dense>
#include <iostream>
#include <vector>

 #include "timer.h"

using namespace Eigen;
using namespace std;

//! \brief Compute the Matrix product $A \times B$ using Strassen's algorithm.
//! \param[in] A Matrix $2^k \times 2^k$
//! \param[in] B Matrix $2^k \times 2^k$
//! \param[out] Matrix product of A and B of dim $2^k \times 2^k$

MatrixXd strassenMatMult(const MatrixXd & A, const MatrixXd & B)
{
    int n=A.rows();
    MatrixXd C(n,n);
    
    if (n==2)
    {
        C<< A(0,0)*B(0,0) + A(0,1)*B(1,0),
            A(0,0)*B(0,1) + A(0,1)*B(1,1),
            A(1,0)*B(0,0) + A(1,1)*B(1,0),
            A(1,0)*B(0,1) + A(1,1)*B(1,1);
        return C;
    }
    
    else
    {   MatrixXd Q0(n/2,n/2),Q1(n/2,n/2),Q2(n/2,n/2),Q3(n/2,n/2),
        Q4(n/2,n/2),Q5(n/2,n/2),Q6(n/2,n/2);
        
        MatrixXd A11=A.topLeftCorner(n/2,n/2);
        MatrixXd A12=A.topRightCorner(n/2,n/2);
        MatrixXd A21=A.bottomLeftCorner(n/2,n/2);
        MatrixXd A22=A.bottomRightCorner(n/2,n/2);
        
        MatrixXd B11=B.topLeftCorner(n/2,n/2);
        MatrixXd B12=B.topRightCorner(n/2,n/2);
        MatrixXd B21=B.bottomLeftCorner(n/2,n/2);
        MatrixXd B22=B.bottomRightCorner(n/2,n/2);
        
        Q0=strassenMatMult(A11+A22,B11+B22);
        Q1=strassenMatMult(A21+A22,B11);
        Q2=strassenMatMult(A11,B12-B22);
        Q3=strassenMatMult(A22,B21-B11);
        Q4=strassenMatMult(A11+A12,B22);
        Q5=strassenMatMult(A21-A11,B11+B12);
        Q6=strassenMatMult(A12-A22,B21+B22);
        
        C<< Q0+Q3-Q4+Q6 ,
        Q2+Q4,
        Q1+Q3,
        Q0+Q2-Q1+Q5;
        return C;
    }
}


int main(void)
{
    srand((unsigned int) time(0));
    
    //check if strassenMatMult works
    int k=2;
    int n=pow(2,k);
    MatrixXd A=MatrixXd::Random(n,n);
    MatrixXd B=MatrixXd::Random(n,n);
    MatrixXd AB(n,n), AxB(n,n);
    AB=strassenMatMult(A,B);
    AxB=A*B;
    cout<<"Using Strassen's method, A*B="<<AB<<endl;
    cout<<"Using standard method, A*B="<<AxB<<endl;
    cout<<"The norm of the error is "<<(AB-AxB).norm()<<endl;
    
    //compare runtimes of strassenMatMult and of direct multiplication
    
    unsigned int repeats = 10;
    timer<> tm_x, tm_strassen;
    std::vector<int> times_x, times_strassen;
    
    for(unsigned int k = 4; k <= 10; k++) {
        tm_x.reset();
        tm_strassen.reset();
        for(unsigned int r = 0; r < repeats; ++r) {
            unsigned int n = pow(2,k);
            A = MatrixXd::Random(n,n);
            B = MatrixXd::Random(n,n);
            MatrixXd AB(n,n);
            
            tm_x.start();
            AB=A*B;
            tm_x.stop();
            
            tm_strassen.start();
            AB=strassenMatMult(A,B);
            tm_strassen.stop();
        }
        std::cout << "The standard matrix multiplication took:       " << tm_x.avg().count() / 1000000. << " ms" << std::endl;
        std::cout << "The Strassen's algorithm took:       " << tm_strassen.avg().count() / 1000000. << " ms" << std::endl;
        
        times_x.push_back( tm_x.avg().count() );
        times_strassen.push_back( tm_strassen.avg().count() );
    }
    
    for(auto it = times_x.begin(); it != times_x.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;
    for(auto it = times_strassen.begin(); it != times_strassen.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;
    
}
