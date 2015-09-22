#include <iostream>
#include <math.h>
#include <Eigen/Dense>

using namespace Eigen

VectorXd arrowmatvec(const VectorXd &d, const VectorXd &a, const VectorXd &x){
  if(d.size() != a.size())
    std::cerr << "Vectors are not the same size!";
  int n = d.size();
  MatrixXd A(n,n);
  A << d.head(n-1).asDiagonal(), a.head(n-1),
       a.head(n-1).transpose(), d.tail(1);
  return A*A*x;
}

int main(){
  
  return 0;
}
