#include <iostream>
#include <math.h>
#include <Eigen/Dense>

using namespace Eigen;

VectorXd arrowmatvec(const VectorXd &d, const VectorXd &a, const VectorXd &x){
  if(d.size() != a.size())
    std::cerr << "Vectors are not the same size!";
  int n = d.size();
  MatrixXd A(n,n);
  // Apparently we can't use the nice << syntax due to asDiagonal
  // screwing it all up.
  A.topLeftCorner(n-1,n-1) = d.head(n-1).asDiagonal();
  A.topRightCorner(1,n-1) = a.head(n-1);
  A.bottomLeftCorner(n-1,1) = a.head(n-1).transpose();
  A.bottomRightCorner(1,1) = d.tail(1);
  return (A*A)*x;
}

VectorXd arrowmatvec2(const VectorXd &d, const VectorXd &a, const VectorXd &x){
  if(d.size() != a.size())
    std::cerr << "Vectors are not the same size!";
  int n = d.size();
  MatrixXd A(n,n);
  A.topLeftCorner(n-1,n-1) = d.head(n-1).asDiagonal();
  A.topRightCorner(1,n-1) = a.head(n-1);
  A.bottomLeftCorner(n-1,1) = a.head(n-1).transpose();
  A.bottomRightCorner(1,1) = d.tail(1);
  return A*(A*x);
}

int main(){
  return 0;
}
