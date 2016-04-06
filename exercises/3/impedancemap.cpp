#include <Eigen/Dense>
#include <Eigen/LU>

#include <iostream>
#include <iomanip>

template <class Matrix>
class ImpedanceMap {
public:
    ImpedanceMap(double R_, double W_) : R(R_), W(W_) {
        // TODO: build and factorize A0 into lu
    };
    
    double operator()(double Rx) {
        // TODO: compute the impedance from voltages and resistance
    };
private:
    Eigen::PartialPivLU< Matrix > lu; //< Store lu decomposition of a for efficiency
    double R, W; //< Resistance R and source voltage W
    Matrix rhs; //< Store rhs vector prescribing sink and source voltages
};

int main(void) {
    ImpedanceMap<Eigen::MatrixXd> IM = ImpedanceMap<Eigen::MatrixXd>(1, 1);
    
    std::cout << std::setw(30) << "Impedance [Ohm]" << std::setw(30) << "R_x [Ohm]" << std::endl;
    for(auto Rx = 1; Rx <= 1024; Rx *= 2) {
        std::cout << std::setw(30) << IM(Rx)        << std::setw(30) << " " << Rx << std::endl;
    }
    
}
