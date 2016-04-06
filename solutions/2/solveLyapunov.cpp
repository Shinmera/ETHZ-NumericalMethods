#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;

//! \brief Compute the matrix C from A.
//! \param[in] A Matrix $n \times n$
//! \param[out] C MatrixXd with $C=A\otimes I+I\otimes A$.

SparseMatrix<double> buildC(const MatrixXd &A)
{
    int n=A.rows();
    
    Eigen::SparseMatrix<double> C(n*n,n*n);
    std::vector<Triplet<double> > triplets;
    MatrixXd I=MatrixXd::Identity(n,n);

    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        if (i==j){
          for (int k1=0; k1<n; k1++) {
            for (int k2=0; k2<n; k2++) {
              Triplet<double>
              triplet(i*n+k1,j*n+k2,A(i,j)*I(k1,k2)+A(k1,k2));
              triplets.push_back(triplet);
             }
          }
        }
        else {
          for (int k=0; k<n ; k++) {
            Triplet<double> triplet(i*n+k,j*n+k,A(i,j));
            triplets.push_back(triplet);
          }
        }
      }
    }
    C.setFromTriplets(triplets.begin(), triplets.end());
    C.makeCompressed();
    return C;
}

//! \brief Solve the Lyapunov system
//! \param[in] A Matrix $n \times n$
//! \param[out] X MatrixXd, the solution.

void solveLyapunov(const MatrixXd &A, MatrixXd &X)
{
    int n=A.rows();
    SparseMatrix <double> C;
    C=buildC(A);
    MatrixXd I=MatrixXd::Identity(n,n);
    VectorXd b(n*n);
    b=Map<MatrixXd>(I.data(),n*n,1);
    VectorXd vecX(n*n);
    SparseLU<SparseMatrix <double> > solver;
    solver.compute(C) ;
    vecX = solver.solve(b);
    X=Map<MatrixXd>(vecX.data(),n,n);
}

int main(){
    
    // test buildC
    int n=5;
    MatrixXd A(n,n);
    A<<10, 2, 3, 4, 5, 6, 20, 8, 9, 1, 1, 2, 30, 4, 5, 6, 7, 8, 20, 0, 1, 2, 3, 4, 10;
    
    SparseMatrix <double> C;
    C=buildC(A);
    cout<<"C= "<<C<<endl;
    
    // solve lynear system
    MatrixXd X(n,n);
    solveLyapunov(A,X);
    cout<<"X= "<<X<<endl;
    MatrixXd I=MatrixXd::Identity(n,n);
    
    // test to verify the solution
    cout<<(A*X+X*A.transpose()-I).norm()<<endl;
}