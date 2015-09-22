#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <chrono>

/*
  Matlab function:

  function Z = houserefl(v)
    n = size(v,1);
    w = v/norm(v);
    u = w + [1;zeros(n-1,1)];
    q = u/norm(u);
    X = eye(n) - 2*q*q';
    Z = X(:,2:end);
  end
 */

using namespace Eigen;

void houserefl(const VectorXd &v, MatrixXd &Z){
  int n = v.cols();
  VectorXd w = v.normalized();
  w(0) += 1;
  w.normalize();
  MatrixXd X = MatrixXd::Identity(n,n) - 2*w*w.transpose();
  Z = X.richtCols(n-1);
}

int main(){
  return 0;
}
