#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

//! \brief Compute F(u), where F is the map associated with the nonlinear circuit
//! \param[in] U VectorXd of size 3
//! \param[in] double alpha, beta parameters
//! \param[in] double Uin, input voltage
//! \param[out] VectorXd, F(U)

VectorXd F(VectorXd U, double alpha, double beta, double Uin)
{
    VectorXd f(3);
    double Ut=.5;
    f << 3*U(0)-U(1)-U(2), 3*U(1)-U(0)-U(2)-Uin, 3*U(2)-U(0)-U(1) + alpha*(exp(beta*(U(2)-Uin)/Ut)-1);
    return f;
}

//! \brief Compute the solution to the nonlinear system by using Newton's iteration
//! \param[in] double alpha, beta parameters
//! \param[in] double Uin, input voltages
//! \param[out] VectorXd Uot, output voltages

void circuit(const double & alpha, const double & beta, const VectorXd & Uin, VectorXd & Uout)
{
    double Ut=.5;
    int n=Uin.size();
    MatrixXd J(3,3);
    VectorXd f(3);
    for (int i=0;i<n;i++){
        srand((unsigned int) time(0));
        VectorXd U=VectorXd::Random(3);
        VectorXd h=VectorXd::Ones(3);
        while (h.cwiseAbs().maxCoeff()>1e-6*U.norm()){
            J << 3, -1, -1, -1, 3, -1, -1, -1, 3+(alpha*beta)/Ut*exp(beta*(U(2)-Uin(i))/Ut);
            f=F(U,alpha,beta,Uin(i));
            h=J.partialPivLu().solve(f);
            U=U-h;
        }
        Uout(i)=U(0);
            // Check correctedness
        cout<<"Error: "<<F(U,alpha,beta,Uin(i)).norm()<<endl;
    }
}


int main(){
    //Test the above function with the input Uin
    int n=20;
    double alpha,beta;
    alpha=8;
    beta=1;
    
    VectorXd Uin=VectorXd::LinSpaced(n,0,20);
    VectorXd Uout(n);
    circuit(alpha,beta,Uin,Uout);
    
    //Display the solutions
    cout<<"The solutions are"<<endl<<Uout<<endl;
    
    // Display the differences of the solutions: the nonlinear effect can be seen from the fact that the vector is not constant.
    VectorXd diff(n-1);
    for (int i=0; i<n-1;i++) {
        diff(i)=Uout(i+1)-Uout(i);
    }
    cout<<"The differences are"<<endl<<diff<<endl;
}
