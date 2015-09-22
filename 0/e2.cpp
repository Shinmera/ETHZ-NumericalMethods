#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <chrono>

using Eigen::MatrixXd;
using namespace std;

MatrixXd strassenMatMult(const MatrixXd &A, const MatrixXd &B){
  int n = A.rows();
  int l = n/2;
  MatrixXd A11 = A.topLeftCorner(l,l);
  MatrixXd A12 = A.topRightCorner(l,l);
  MatrixXd A21 = A.bottomLeftCorner(l,l);
  MatrixXd A22 = A.bottomRightCorner(l,l);
  MatrixXd B11 = B.topLeftCorner(l,l);
  MatrixXd B12 = B.topRightCorner(l,l);
  MatrixXd B21 = B.bottomLeftCorner(l,l);
  MatrixXd B22 = B.bottomRightCorner(l,l);

  MatrixXd Q0(l,l),Q1(l,l),Q2(l,l),Q3(l,l),Q4(l,l),Q5(l,l),Q6(l,l);

  if(l==1){
    Q0 = (A11+A22)*(B11+B22);
    Q1 = (A21+A22)*B11;
    Q2 = A11*(B12-B22);
    Q3 = A22*(-B11+B21);
    Q4 = (A11+A12)*B22;
    Q5 = (-A11+A21)*(B11+B12);
    Q6 = (A12-A22)*(B21+B22);
  }else{
    Q0 = strassenMatMult((A11+A22),(B11+B22));
    Q1 = strassenMatMult((A21+A22),B11);
    Q2 = strassenMatMult(A11,(B12-B22));
    Q3 = strassenMatMult(A22,(-B11+B21));
    Q4 = strassenMatMult((A11+A12),B22);
    Q5 = strassenMatMult((-A11+A21),(B11+B12));
    Q6 = strassenMatMult((A12-A22),(B21+B22));
  }

  MatrixXd C(n,n);
  C.topLeftCorner(l,l) = Q0+Q3-Q4+Q6;
  C.topRightCorner(l,l) = Q2+Q4;
  C.bottomLeftCorner(l,l) = Q1+Q4;
  C.bottomRightCorner(l,l) = Q0-Q1+Q2+Q5;
  return C;
}

void e2b(){
  MatrixXd A = MatrixXd::Random(4,4);
  MatrixXd B = MatrixXd::Random(4,4);
  MatrixXd ABS = strassenMatMult(A,B);
  MatrixXd ABE = A*B;

  cout << "Strassen: " << endl;
  cout << ABS << endl;
  cout << "Error: " << endl;
  cout << (ABS-ABE) << endl;
}

void e2c(){
  for(int k = 2; k<= 8; ++k){
    unsigned int n = pow(2,k);
    cout << "N: " << n << endl;

    MatrixXd A = MatrixXd::Random(n,n);
    MatrixXd B = MatrixXd::Random(n,n);
    
    auto ts1 = chrono::high_resolution_clock::now();
    MatrixXd ABS = strassenMatMult(A,B);
    auto ts2 = chrono::high_resolution_clock::now();
    auto te1 = chrono::high_resolution_clock::now();
    MatrixXd ABE = A*B;
    auto te2 = chrono::high_resolution_clock::now();
    cout << "Strassen: " << chrono::duration_cast<chrono::milliseconds>(ts2-ts1).count() << " ms" << endl;
    cout << "Eigen:    " << chrono::duration_cast<chrono::milliseconds>(te2-te1).count() << " ms" << endl;
    cout << endl;
  }
}

int main(void) {
  e2b();
  e2c();
  return 0;
}
