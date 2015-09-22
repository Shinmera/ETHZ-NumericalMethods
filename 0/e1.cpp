#include <iostream>
#include <Eigen/Dense>

/*
Original Matlab Code 1.5.3:

  function Q = gramschmidt(A)
  % Gram-Schmidt orthogonalization of column vectors
  % Arguments: Matrix A passes vectors in its columns
  % Return values: orthornormal system in columns of matrix Q
  [n,k] = size(A);                               % Get number k of vectors and dimension n of space
  Q = A(:,1)/norm(A(:,1));                       % First basis vector
  for j=2:k
      q = A(:,j) - Q*(Q’*A(:,j));                % Orthogonal projection; loop-free implementation
      nq = norm(q);                              % Check premature termination
      if (nq < (1E-9) *norm(A(:,j))), break; end % Safe check for == 0
      Q = [Q,q/nq];                              % Add new basis vector as another column of Q
  end
*/

template <class Matrix>
Matrix gramschmidt(const Matrix &A){
  Matrix Q = A;
  Q.col(0).normalize();
  
  for (int j=1; j< A.cols(); ++j){
    Q.col(j) -= Q.leftCols(j) * (Q.leftCols(j).transpose() * A.col(j));

    if(Q.col(j).norm() < (10e-14 * A.col(j).norm())){
      std::cerr << "Can't normalize!" << std::endl;
    }

    Q.col(j).normalize();
  }
  return Q;
}


int main(void) {
  // Ortho test
  unsigned int n = 9;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
  Eigen::MatrixXd Q = gramschmidt(A);
    
  // Output should be idenity matrix
  std::cout << Q*Q.transpose() << std::endl;
}
